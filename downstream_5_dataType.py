"""
Downstream 5: Converts the PAUP files (output of DS3), 
cleans up the data types, and applies coverage filtering

Inputs: 
    1) paup file of choice (output of DS3)
    2) host information in the form of a CSC
    3) name of the new file without extension
    4) number above which coverage if chosen 
    5) coverage file (output of DS4)
"""

import pandas as pd 
import sys
import numpy as np
from numpy.random import choice

def importCoverageDict(filename, cleanLabels):
    coverageDict = {}
    with open(filename, 'r') as f:
        for line in f:
            temp = {}
            splitLine = line.split(";")
            host = splitLine[0].replace(" ", "_")
            data = splitLine[1:]
            for d in data:
                if ":" in d:
                    splitD = d.split(":")
                    bact = splitD[0].replace(" ", "_").replace("-", "")
                    if bact in cleanLabels:
                        bIndex = cleanLabels.index(bact)
                        temp[bIndex] = int(splitD[1])
                    #else:
                        #print(host, bact)
            coverageDict[host] = temp
    return coverageDict

def writingNexus(d, dataType, chars, location):
    filename = "paupFiles/" + location + "_" + dataType + ".nex"
    print(len(d))
    with open(filename, 'w') as f:
        f.write("#NEXUS\n")
        f.write("begin data;\n")
        f.write("dimensions ntax=%i nchar=%i;\n" %(len(d), len(chars)))
        charString = "charlabels " + " ".join(chars) + ";\n"
        f.write(charString)
        f.write("matrix\n")
        for key in d:
            f.write(key + " " + " ".join(str(i) for i in d[key]) + "\n")
        f.write(";\n")
        f.write("end;\n")

def writingHosts(d, dataType, location):
    filename = "hostFiles/" + location + "_" + dataType + ".txt"
    with open(filename, 'w') as f:
        for key in d:
            temp = key.replace("_", " ")
            f.write(temp + '\n')

def mergeNames(original, add): 
    for i in range(len(original)):
        if add[i] == 1 and original[i] == 0:
            original[i] == 1
    return original

def writingFiles(hostFile, newDict,cleanLabels, location, threshold, coverageDict):
    hostsdf = pd.read_csv(hostFile)
    # Mosquito: Anopheles_sp
    updatedDict = {}

    for key in newDict:
        if key in hostsdf['Newick_label'].values: 
            tempIndex = hostsdf[hostsdf['Newick_label']==key].index.tolist()
            tempDF = hostsdf.iloc[tempIndex[0]]
            if tempDF['correctedName'] == np.nan:
                if key not in updatedDict:
                    updatedDict[key] = newDict[key]
                else:
                    updatedDict[key] = mergeNames(updatedDict[key], newDict[key])
            else:
                if tempDF['correctedName'] in newDict.keys():
                    updatedDict[tempDF['correctedName']] = mergeNames(newDict[tempDF['correctedName']], newDict[key])
                    print("merged this key %s with %s" %(key, tempDF['correctedName'])) 
                else:
                    updatedDict[key] = newDict[key]
        else:
            print("removing this key", key)
    
    num = threshold
    print()
    print("Coverage updates for the following:")
    for key in updatedDict:
        if sum(updatedDict[key]) > num:
            print(key)
            print(coverageDict[key])
            indices = []
            probs = []
            totalProb = 0
            for i in range(len(updatedDict[key])):
                if updatedDict[key][i] == 1:
                    #indices.append(i)
                    if i in coverageDict[key]:
                        indices.append(i)
                        totalProb += coverageDict[key][i]
                    else:
                        print(i, cleanLabels[i])
            for index in indices:
                probs.append( float(coverageDict[key][index]) / float(totalProb) )           
            if len(indices) < num:
                newIndices = indices    
            else:
                newIndices = choice(indices, size=num, p=probs, replace=False)
            print(key, len(indices), newIndices)
            for n in newIndices:
                print(cleanLabels[n])
            print()
            sampled = []
            for j in range(len(updatedDict[key])):
                if j in newIndices:
                    sampled.append(1)
                else:
                    sampled.append(0)
            updatedDict[key] = sampled
    
    df = pd.DataFrame(updatedDict).T
    df.columns=cleanLabels
    df.loc["Total"] = df.sum()
    cols = []
    for i in range(len(df.loc["Total"])):
        val = df.loc["Total"][i]
        if val == 0:
            cols.append(df.columns[i])
    newDF = df.drop(cols,1)
    print(newDF.shape)
    print()
    finalLabels = newDF.columns
    d = {}
    for i in range(len(newDF)):
        if newDF.iloc[i].name != "Total":
            d[newDF.iloc[i].name] = newDF.iloc[i].values
    writingNexus(d, "updated", finalLabels,location)
    writingHosts(d, "updated", location)
    return

def openPaupFile(filename):
    with open(filename, 'r') as f: 
        lines = []
        m = False
        for line in f: 
            if 'charlabels' in line: 
                labels = line
            elif 'matrix' in line or m:
                m = True
                lines.append(line.rstrip())
            elif line == ';\n': 
                m = False
    lines = lines[1:-2]
    print("Length of Lines: ", len(lines))

    splitLabels = labels.split(" ")
    splitLabels[-1] = splitLabels[-1].rstrip().replace(';','')
    cleanLabels = splitLabels[1:]

    lines = [x.rstrip() for x in lines]
    linesDict = {}
    for val in lines: 
        splitVal = val.split(" ")
        if ';' in splitVal or '#nexus' in splitVal or 'begin' in splitVal or 'dimensions' in splitVal or 'matrix' in splitVal or 'end;' in splitVal: 
            continue
        else: 
            linesDict[splitVal[0]] = [int(y) for y in splitVal[1:]]

    return linesDict, cleanLabels

def main():
    linesDict, labels = openPaupFile(sys.argv[1])
    coverageDict = importCoverageDict(sys.argv[5], labels)
    writingFiles(sys.argv[2], linesDict,labels, sys.argv[3], int(sys.argv[4]), coverageDict)

main()
