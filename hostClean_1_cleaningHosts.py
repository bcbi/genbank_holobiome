"""
Cleaning Hosts
ONLY if you are updating cleaned hosts annotations
"""

import pandas as pd
import sys
sys.path.append("/users/vramanan/.local/bin")

removeWords =["ISOLATE", "ADULT", "Adult", "Spot", "spot", "View", "'s", "SP", "TFF2", \
        "wt", "Allele", "protein", ", human", "extract", "dietary", "phosphoethanolamine", \
       "palmidrol", "APACHE", "Apache", "apache", "-", "LI", "AR", "AKR1B1", "Li+", \
       "Androgen", "Receptor", "Positive", "Plant", "Cone", "Scaling", "Gut", "GUT", \
       "Lung", "LUNG", "Strain", "STRAIN", "root", "ROOT", "Root", "heirarchal", "nodule",\
       "NODULE", "Nodule", "body", "part", "Tooth", "structure", "Roots", "Type", "type",\
       "strain", "TYPE", "BAHAMAS", "Bahamas", "FALSE", "false", "False", "Common", \
       "common", "subsp.", "INTESTINE", "Intestines", "Measurement", "Third", "stage", \
       "LARVA", "Larva", "larva", "Foregut", "Foregu", "foregut", "Primitive", "Larvae",\
       "Laboratory", "cultur", "Group", "Two", "MARK", "mark", "Mark", "clone", "genotype",\
        "Pool", "pool", "Environmental", "action", "Sample", "Stage", "stage", "Phase", \
        "Tumor", "BREED", "breed", "group", "Breed", "NOS", "NRF", "MALE", "FEMALE", "unknown",\
        "UNKNOWN"]

def removeParanthesis(value): 
    """
    Removes paranthesis from a given value
    """
    start = value.find('(')
    end = value.find(")")
    new = value[0:start] + value[end:len(value)-1]
    return new

def removeDigits(s):
    """
    Removes all numerical values from given string s
    """
    answer = []
    for char in s:
        if (not char.isdigit()):
            answer.append(char)
    return ''.join(answer)

def speciesDictLoad(): 
    """
    Loads in speciesNames.txt
    """
    with open("/users/vramanan/genbank_database_pipeline/speciesNames.txt", 'r') as f: 
        speciesDict = {}
        for line in f: 
            splitLine = line.split(":")
            scName = splitLine[0]
            values = splitLine[1].split(',')
            valuesStripped = [val.rstrip() for val in values]
            speciesDict[scName] = valuesStripped
    return speciesDict

def oneWordUpdates(cleanHost, speciesDict, data):
    """
    Updates all oneWord potential values into species names from speciesNames.txt
    """

    oneWord = []
    for d in cleanHost: 
        if len(d.split()) == 1 and d.strip() not in oneWord: 
            oneWord.append(d.strip())

    scNameDict = {}
    for val in oneWord:
        for scName in speciesDict: 
            if any(word == val for word in speciesDict[scName]):
                scNameDict[val] = scName

    newHost = []
    for i in range(len(cleanHost)):
        val = cleanHost[i].strip()
        if val in scNameDict.keys(): 
            newHost.append(scNameDict[val])
        elif val == "" or val == " ":
            oldValue = data['cleanScientificName'][i]
            # Fix: “Horse”, “Canine”, “Chicken”, “Human”, “Cow”, “Mouse”, “Cattle"
            if "horse" in oldValue.lower():
                newHost.append("Equus caballus")
            elif "canine" in oldValue.lower():
                newHost.append("Canis familiar")
            elif "chicken" in oldValue.lower():
                newHost.append("Gallus gallus")
            elif "human" in oldValue.lower() or "homo " in oldValue.lower():
                newHost.append("Homo sapiens")
            elif "cow" in oldValue.lower() or "cattle" in oldValue.lower():
                newHost.append("Bos taurus")
            elif "mouse" in oldValue.lower() or oldValue=="Mouse," or oldValue=="ob/ob Mouse":
                newHost.append("Mus musculus")
            elif "goat" in oldValue.lower(): 
                newHost.append("Capra aegagrus hircus")
            elif oldValue=="Fly" or oldValue=="Fly,":
                newHost.append("Drosophila melanogaster")
            elif oldValue=="Rat":
                newHost.append("Rattus")
            elif "porcine" in oldValue.lower():
                newHost.append("Sus scrofa")
            else:
                newHost.append(oldValue)
        else:
            newHost.append(val)

    return newHost

def cleaning(data, taxFile): 
    """
    Executes whole cleaning on hosts
    """

    cleanHost = []
    total = 0
    speciesDict = speciesDictLoad()
    data.columns = ["locus", "originalName", "cleanScientificName", "species"]
    totalLength = float(len(data['cleanScientificName']))

    for value in data['cleanScientificName']:
        # remove unnecessary words
        if any(word in value for word in removeWords):
            splitVal = value.split()
            new = ' '.join([word for word in splitVal if word not in removeWords])
            newValue = new.strip()
        else:
            newValue = value.strip()
        
        # check for scientific names
        for scName in speciesDict.keys():
            if scName in newValue: 
                newValue = scName
                break
            for val in speciesDict[scName]:
                if val in newValue: 
                    newValue = scName
                    break
        
        # updates for () or other symbols 
        if "(" in newValue: 
            new = removeParanthesis(newValue)
            if '(' in new or ')' in new: 
                new = new.replace('(', '').replace(')','')
            new = removeDigits(new)
            new = new.replace("%", '')
            new = new.replace("*", '')
            new = new.replace("<", '')
            new = new.replace(">", '')
            newValue = new
        elif "," in newValue: 
            new = newValue.replace(",", '')
            newValue = new
        
        # Updating capitalization
        if newValue.isupper(): 
            caps = newValue.lower()
            new = caps.capitalize()
            newValue = new
        elif newValue.islower():
            new = newValue.capitalize()
            newValue = new
            
        # Manual updates
        if "rice" in newValue or "Oryza" in newValue or "Rice" in newValue: 
            newVal = "Oryza sativa"
        elif "Broiler" in newValue:
            newVal = "Chicken"
        elif "Magnolia grandiflora" in newValue:
            newVal = "Magnolia grandiflora"
        elif "Drosophila melanogaster" in newValue or newValue=="Fly" or newValue=="Fly,":
            newVal = "Drosophila melanogaster"
        elif "Mus musculus" in newValue or "Musculus" in newValue or newValue == "Mouse" or newValue=="Mouse," or newValue=="ob/ob Mouse":
            newVal = "Mus musculus"
        elif "Seal" in newValue:
            newVal = "Seal"
        elif "Lutzomyia" in newValue:
            newVal = "Lutzomyia"
        elif "Aphid" in newValue or "APHID" in newValue or "aphid" in newValue:
            newVal = "Aphidoidea"
        elif "Weevil" in newValue or "WEEVIL" in newValue or "weevil" in newValue:
            newVal = "Curculionoidea"
        elif "Cervus nippon" in newValue: 
            newVal = "Cervus nippon"
        elif "Vitis vinifera" in newValue:
            newVal = "Vitis vinifera"
        elif "Muskoxen" in newValue:
            newVal = "Muskoxen ovibus"
        elif "CHICKEN" in newValue:
            newVal = "Gallus gallus"
        elif "Nematoda" in newValue or "Nematod" in newValue or "nematoda" in newValue or "nematode" in newValue:
            newVal = "Nematoda"
        elif "Wallaby" in newValue or "Wallabies" in newValue:
            newVal = "Wallaby"
        elif "BUFFALO" in newValue:
            newVal = "Buffalo"
        elif newValue=="Ape" or newValue=="Apes":
            newVal = "Primate"
        elif "Panulirus ornatus" in newValue:
            newValue = "Panulirus ornatus"
        elif "Rattus" in newValue or newValue == "Rat" or newValue == "rat":
            if "norvegicus" in newValue:
                newVal = "Rattus norvegicus"
            elif "Rattus Rattus" in newValue:
                newVal = "Rattus rattus"
            elif newValue == "Rat" or "Rattus sp." in newValue:
                newVal = "Rattus"
            else:
                newVal = newValue
        elif "Homo sapiens" in newValue or "homo " in newValue.lower():
            newVal = "Homo sapiens"
        elif "porcine" in newValue.lower():
            newVal = "Sus scrofa"
        elif "capra hircus" in newValue.lower():
            newVal = "Capra aegagrus hircus"
        else:
            newVal = newValue
        
        cleanHost.append(newVal)
        total += 1

        if total % 1000000 == 0: 
            print("Percent Completed: %0.2f" %(float(total)/totalLength), flush=True)

    newHost = oneWordUpdates(cleanHost, speciesDict, data)
    newHost, phylum, classCol, order, family = removeBacteria(newHost, taxFile)

    data['updatedScientificName'] = newHost
    data['phylum'] = phylum
    data['class'] = classCol
    data['order'] = order
    data['family'] = family
    return data

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

def removeBacteria(newHost, taxFile):
    """
    Removes bacterial hosts and makes them empty string
    """
    bacteria = ["streptococcus", "bacillus", "bacteria", "bacterium", "bacterio", "firmicutes",\
                "actinobacteria", "proteobacteria", "bacteroidetes", "escherichia", "pseudomonas", \
                "klebsiella", "bifidobacterium", "lactobacillus"]
    taxonomyDict = {}
    with open(taxFile, 'r') as f:
        for line in f:
            split = line.rstrip().split("\t")
            key = split[0]
            vals = split[1:]
            keyDict = {}
            for value in vals:
                splitVal = value.split(":")
                keyDict[splitVal[0]] = splitVal[1]
            taxonomyDict[key] = keyDict
    print(len(taxonomyDict), flush=True)
    print("removing bacteria", flush=True)
    add = {'phylum':[],'class':[], 'order':[], 'family':[], 'genus':[], 'resolution':[]}
    for i in range(len(newHost)):
        val = newHost[i]
        tax = taxonomyDict[val]
        add = assigntoDict(add, tax) 
        if any(bact in add['phylum'][i] for bact in bacteria):
            print(val, add['phylum'][i], flush=True)
            newHost[i] = "Bacteria"

    return newHost, add['phylum'], add['class'], add['order'], add['family']

def main():

    data = pd.read_csv(sys.argv[1])
    newData = cleaning(data, sys.argv[2])
    newData.to_csv("hostData/updated_cleanedhost.csv", index=False)

main()
