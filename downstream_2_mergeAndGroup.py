"""
midstream 2
"""
import pandas as pd
import sys
from Bio import Entrez
sys.path.append("/users/vramanan/.local/bin")

Entrez.email = "vivek_ramanan@brown.edu" 
Entrez.sleep_between_tries = 15
Entrez.api_key = "f58d32928f53b067d24ec5ddee6c7a503b08"

def EntrezSearch(searchTerm):
    handle = Entrez.esearch(db="taxonomy",term=searchTerm,retmode = "xml")
    record = Entrez.read(handle)
    return record

def EndrezFetch(record, term):

    speciesDict = {'TM7':{'ScientificName':'Candidatus Saccharibacteria','Rank':'phylum'},
               'Measles virus':{'TaxId': '11229', 'ScientificName': 'Morbillivirus', 'Rank': 'genus'}}

    try: 
        searchID = record['IdList'][0]
    except:
        if "PhraseNotFound" in record['ErrorList']:
            for spec in speciesDict: 
                if spec in term: 
                    return [speciesDict[spec]]
            return [{'FAIL':'Phrase Not Found','ScientificName':term}]
    search = Entrez.efetch(id = searchID, db = "taxonomy", retmode = "xml")
    output = Entrez.read(search)
    try:
        d = output[0]['LineageEx']
    except:
        return [{'FAIL':'Lineage Ex Error','ScientificName':term}]
    for i in range(len(d)-1,0,-1):
        val = d[i]
        if val['Rank'] == 'genus':
            return [val, d[i-1]]
        elif val['Rank'] == 'family':
            return [val]
        elif val['Rank'] == 'order':
            return [val]
        elif val['Rank'] == 'phylum':
            return [val]
    
    return [{'FAIL':'Not enough resolution','ScientificName':term}]

def updateSpecies(df): 

    speciesClean = []
    for val in df['species']:
        valSplit = val.split("(")[0].strip()
        valSplit = valSplit.replace("[",'').replace("]",'')
        speciesClean.append(valSplit)
    speciesCleanSeries = pd.Series(speciesClean)

    taxonomyDict = {}
    for d in speciesCleanSeries.value_counts().index:
        record = EntrezSearch(d)
        output = EndrezFetch(record, d)
        taxonomyDict[d] = output

    genus = []
    family = []
    resolution = []
    for val in speciesCleanSeries:
        tax = taxonomyDict[val]
        if len(tax) == 2:
            genus.append(tax[0]['ScientificName'])
            family.append(tax[1]['ScientificName'])
            resolution.append('genus')
        else: 
            # includes both fails and single values (family, order, phyla)
            genus.append(tax[0]['ScientificName'])
            family.append(tax[0]['ScientificName'])
            if 'FAIL' in tax[0]:
                resolution.append('Fail')
            else:
                resolution.append(tax[0]['Rank'])
    
    df['speciesGenus'] = genus
    df['speciesFamily'] = family

    return df

def updateGroupings(df, groupingFile): 

    groupingsDict = {}
    with open(groupingFile, 'r') as f: 
        for line in f:
            splitLine = line.split(":")
            key = splitLine[0]
            valueSplit = splitLine[1].split(",")
            values = [val.rstrip() for val in valueSplit]
            groupingsDict[key] = values

    print(groupingsDict)

    groups = []
    for val in df['preferred']: 
        tempGrouping = []
        for group in groupingsDict: 
            if any(word in val for word in groupingsDict[group]): 
                tempGrouping.append(group)
        if len(tempGrouping) >= 1: 
            groups.append(tempGrouping[0])
        else: 
            groups.append("MISC")

    print(len(groups))
    print(len(df['preferred']))
                
    df['group'] = groups

    return df

def combineDatabases(file1, file2, groupingFile): 
    df1 = pd.read_csv(file1, low_memory=False)
    df2 = pd.read_csv(file2, low_memory=False)
    df1['DataType'] = 'isolation'
    df2['DataType'] = 'tissue'

    merged = pd.concat([df1, df2])
    merged_Groups = updateGroupings(merged, groupingFile)
    #merged_Species = updateSpecies(merged_Groups)

    return merged_Groups

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    groupingFile = sys.argv[3]
    merged = combineDatabases(file1, file2, groupingFile)
    merged.to_csv("merged_Host_Bacteria.csv",index=False)


main()
