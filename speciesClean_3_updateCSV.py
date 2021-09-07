"""
Updating Species to Genus/Family
Vivek Ramanan
"""

import pandas as pd
import sys, time
from Bio import Entrez
sys.path.append("/home/vramanan/.local/bin")

def assigntoDict(add, tax):
    """
    Function: assign to Dict - takes a singular value from taxonomyDict to parse into correct format
    @param add: parsed values go into this dict
    @param tax: the particular value to parse
    """
    if tax['resolution'] == 'fail':
        for key in add:
            add[key].append(tax['FailTerm'])
    elif len(tax) == 6: 
        for key in tax: 
            add[key].append(tax[key])
    else: 
        for val in add:
            if val in tax:
                add[val].append(tax[val])
            else:
                resolution = tax['resolution']
                add[val].append(tax[resolution])
    return add

def updateSpecies(df, speciesFile): 
    """
    Function: updating Species
    @param df: dataframe that is input, has column 'species' 
    @return updated cleaned species df
    """

    print("Starting cleaning of species.", flush=True)
    taxonomyDict = {}
    with open(speciesFile, 'r') as f:
        for line in f:
            split = line.rstrip().split("\t")
            key = split[0]
            vals = split[1:]
            if "Sendai virus" in key:
                print(key,flush=True)
                vals = ["genus:Respirovirus","family:Paramyxoviridae","order:Mononegavirales",\
                        "class:Monjiviricetes","phylum:Negarnaviricota","resolution:genus"]
            elif "Bastrovirus" in key:
                print(key,flush=True)
                vals = ["family:Astroviridae","order:Stellavirales","class:Stelpaviricetes",\
                        "phylum:Pisuviricota","resolution:family"]
            keyDict = {}
            for value in vals:
                splitVal = value.split(":")
                keyDict[splitVal[0]] = splitVal[1]
            taxonomyDict[key] = keyDict
    print(len(taxonomyDict))
    
    add = {'genus':[],'family':[],'order':[],'class':[],'phylum':[],'resolution':[]}
    for val in df['species']:
        tax = taxonomyDict[val]
        add = assigntoDict(add, tax)

    print("Assigned to dictionary", flush=True)
    
    df['speciesGenus'] = add['genus']
    df['speciesFamily'] = add['family']
    df['speciesOrder'] = add['order']
    df['speciesClass'] = add['class']
    df['speciesPhylum'] = add['phylum']
    df['speciesResolution'] = add['resolution']
    return df

def writeDict(f, val, output):
    """
    Writes the taxonomy dictionary to file so it's saved
    @param f: taxonomyDict.txt
    @param val: the original search term
    @param output: dictionary from the Entrez output
    """
    stringToWrite = val + '\t'
    for k in output:
        temp = k + ":" + output[k] + '\t'
        stringToWrite += temp
    f.write(stringToWrite+'\n')
    return

def startSpecies(dfFile, sFile): 
    """
    Starts the process - reads in file and runs updateSpecies
    """
    print("Reading in file", flush=True)
    df = pd.read_csv(dfFile,low_memory=False)
    merged_Species = updateSpecies(df, sFile)

    return merged_Species

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    merged = startSpecies(file1, file2)
    merged.to_csv("merged_Host_Bacteria_Species.csv",index=False)


main()
