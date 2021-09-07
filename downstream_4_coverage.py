"""
Downstream 4: Creates a coverage file for DS5 to use 
Input: 1) hostFiles/<group>_coverage.txt (one output of DS3)
        2) name of new coverage file to use in DS5
Author: Vivek Ramanan

"""

import sys, string

def turnListDict(line):
    table = str.maketrans(dict.fromkeys(string.punctuation))  # OR {key: None for key in string.punctuation}
    
    split = line.split(",")
    retDict = {}
    for s in split:
        val = s.translate(table)
        val = val.strip()
        if val in retDict:
            retDict[val] += 1
        else:
            retDict[val] = 1
    return retDict

def revDict(microDict):
    reverseDict = {}
    for key in microDict:
        for val in microDict[key]:
            if val in reverseDict:
                reverseDict[val][key] = microDict[key][val]
            else:
                reverseDict[val] = {key:microDict[key][val]}
    return reverseDict

def writeRevDict(r, filename):
    with open(filename, 'w') as f:
        for key in r: 
            line = key + ";"
            for val in r[key]:
                line += val + ":" + str(r[key][val]) + ";"
            f.write(line + '\n')

def main():
    filename = sys.argv[2]
    with open(sys.argv[1], 'r') as f:
        micro = ""
        microDict = {}
        for line in f:
            if line.startswith("["):
                microDict[micro] = turnListDict(line)
            else:
                micro = line.rstrip()
    r = revDict(microDict)
    writeRevDict(r, filename)


main()
