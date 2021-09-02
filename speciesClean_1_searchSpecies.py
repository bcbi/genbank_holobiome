"""
Updating Species to Genus/Family
Vivek Ramanan
"""

import pandas as pd
import sys, time
from Bio import Entrez
sys.path.append("/home/vramanan/.local/bin")

Entrez.email = "vivek_ramanan@brown.edu" 
Entrez.sleep_between_tries = 30
Entrez.api_key = "f58d32928f53b067d24ec5ddee6c7a503b08"

def EntrezSearch(searchTerm):
    """
    Runs initial search through Entrez for @param searchTerm
    @return record: object of Entrez search
    """
    try: 
        handle = Entrez.esearch(db="taxonomy",term=searchTerm,retmode = "xml")
    except:
        time.sleep(5)
        handle = Entrez.esearch(db="taxonomy",term=searchTerm,retmode = "xml")
    try: 
        record = Entrez.read(handle)
    except:
        time.sleep(5)
        handle = Entrez.esearch(db="taxonomy",term=searchTerm,retmode = "xml")
        try: 
            record = Entrez.read(handle)
        except:
            time.sleep(5)
            handle = Entrez.esearch(db="taxonomy",term=searchTerm,retmode="xml")
            record = Entrez.read(handle)
    
    return record

def EntrezFetch(record, term):
    """
    Runs the Entrez Fetch process
    @param record: output from Entrez Search function
    @param term: search term
    @return termDict from cleanOutput
    """

    speciesDict = {'TM7':{'phylum':'Candidatus Saccharibacteria','resolution':'phylum'},\
        'Measles virus':{'genus': 'Morbillivirus', 'family':'Paramyxoviridae','order': 'Mononegavirales',\
                          'class': 'Monjiviricetes','phylum': 'Negarnaviricota','resolution': 'genus'},
        'Influenza A': {'genus': 'Alphainfluenzavirus', 'family':'Orthomyxoviridae','order': 'Articulavirales',\
                          'class': 'Insthoviricetes','phylum': 'Negarnaviricota','resolution': 'genus'},
        'Bacillus megaterium':{'genus':'Bacillus', 'family':'Bacillaceae', 'order':'Bacillales', \
                          'class':'Bacilli', 'phylum':'Firmicutes', 'resolution':'genus'},
        'Lactobacillus zeae':{'genus':'Lactobacillus','family':'Lactobacillaceae','order':'Lactobacillales',\
                          'class':'Bacilli','phylum':'Firmicutes','resolution':'genus'}
              }

    try: 
        searchID = record['IdList'][0]
    except:
        if "PhraseNotFound" in record['ErrorList']:
            for spec in speciesDict: 
                if spec in term: 
                    return speciesDict[spec]
            return {'resolution':'fail','FailTerm':term}
    try: 
        search = Entrez.efetch(id = searchID, db = "taxonomy", retmode = "xml")
    except:
        time.sleep(5)
        search = Entrez.efetch(id = searchID, db = "taxonomy", retmode = "xml")
    try: 
        output = Entrez.read(search)
    except:
        time.sleep(5)
        output = Entrez.read(search)
    try:
        d = output[0]['LineageEx']
    except:
        return {'resolution':'fail','FailTerm':term}
    termDict = cleanOutput(d)
    if termDict['resolution'] != 'fail': 
        return termDict
    else: 
        return {'resolution':'fail','FailTerm':term}

def cleanOutput(d):
    """
    Additionaly cleaning step for outputs of the Entrez fetch section
    @param d: output of Entrez.fetch
    @return tempDict
    """
    taxRanks = ['phylum','class','order','family','genus']
    tempDict = {}
    resolutionOptions = []
    for val in d:
        try:  
            if val['Rank'] in taxRanks:
                tempDict[val['Rank']] = val['ScientificName']
                resolutionOptions.append(val['Rank'])
        except:
            continue
    if len(resolutionOptions) != 0:
        tempDict['resolution'] = resolutionOptions[len(resolutionOptions) - 1]
    else:
        tempDict['resolution'] = 'fail' 
    return tempDict

def assigntoDict(add, tax):
    """
    Function: assign to Dict - takes a singular value from taxonomyDict to 
    parse into correct format
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

def updateSpecies(s): 
    """
    Function: updating Species
    """

    print("Length of Clean Species List: %i" %(len(s)), flush=True)

    taxFile = "taxonomyDict.txt"
    taxonomyDict = {}
    total = 0
    f = open(taxFile, 'a')
    totalLength = float(len(s))
    for val in s:
        valSplit = val.split("(")[0].strip()
        d = valSplit.replace("[",'').replace("]",'').replace("/",'')
        record = EntrezSearch(d)
        output = EntrezFetch(record, d)
        if output['resolution'] == 'fail':
            sp = d.split(" ")[0:2]
            d = " ".join(sp)
            if "HIV-1" in d or "virus" in sp[0]:
                d = sp[0]
            elif "kobuvirus" in d.lower():
                d = "Kobuvirus"
            elif "bastrovirus" in d.lower(): 
                d = "Bastrovirus"
            elif "Caliciviridae" in d:
                d = "Caliciviridae"
            elif "tailam virus" in d.lower():
                d = "Tailam virus"
            elif "escherichiashigella" in d.lower():
                d = "Escherichia coli"
            record = EntrezSearch(d)
            output = EntrezFetch(record, d)
            if output['resolution'] == 'fail':
                if "bacterium" in d.lower():
                    output['FailTerm'] = "unclassified bacterium"
                elif "Virus PK" in d or "Virus NG" in d or "Virus Chimp" in d or "Phage" in d:
                    output['FailTerm'] = "unclassified virus"
                elif "uncultured" in d or "unidentified" in d:
                    value = output['FailTerm']
                    output['FailTerm'] = value.replace("unidentified","unclassified").replace("uncultured","unclassified") 
        taxonomyDict[val] = output
        writeDict(f, val, output)
        f.flush()
        total += 1
        if total % 1000 == 0:
            print("%0.2f done" %(float(total)/totalLength), flush=True) 
    f.close()
    print("Length of Taxonomy Dictionary: ",len(taxonomyDict),flush=True)
    return

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

def startSpecies(dfFile): 
    """
    Starts the process - reads in file and runs updateSpecies
    """
    print("Reading in file", flush=True)
    speciesNames = []
    with open(dfFile, 'r') as f:
        for line in f: 
            speciesNames.append(line.rstrip())
    updateSpecies(speciesNames)

    return 

def main():
    startSpecies(sys.argv[1])

main()
