"""
PAUP Processing Functions
"""

def mergeData(l1, l2):
    """
    Merges two lists of presence absence data
    @param l1, l2 lists of 0/1 int data
    @return merged list
    """
    newList = []
    for i in range(len(l1)):
        if l1[i] + l2[i] < 1: 
            newList.append(0)
        else:
            newList.append(1)
    return newList

def updateData(betaData, hostUpdateFile):
    """
    Runs the updating of hosts from hostUpdateFile for naming conventions 
    and merging those with same corrected name
    ** file must have the following columns: Newick_label, correctedName
    @param betaData: dictionary from PAUPprocess, uncleaned
    @param hostUpdateFile: FULL pathway for hostUpdate File
    @return newBetaData: cleaned dictionary
    """
    hostDF2 = pd.read_csv(hostUpdateFile)
    newBetaData = {}
    for b in betaData:
        if b in hostDF2['Newick_label'].values:
            correctedName = hostDF2.loc[hostDF2['Newick_label'] == b]['correctedName'].values[0]
            if type(correctedName)== str: 
                if correctedName in newBetaData: 
                    newBetaData[correctedName] = mergeData(newBetaData[correctedName],betaData[b])
                else: 
                    newBetaData[correctedName] = betaData[b]
            else:
                newBetaData[b] = betaData[b]
    return newBetaData

def PAUPprocess(filename, hostUpdateFile): 
    """
    Processes PAUP files into dictionary while also merging data based on hostUpdates file
    @param filename: FULL pathway to paup file
    @param hostUpdateFile: FULL pathway to host update CSV file
    @return newLinesDict (dictionary where keys:hosts and values:0/1 data arrays), cleanLabels (column labels)
    """
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
    print("Original Length of Lines: ", len(lines))
    
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
    newLinesDict = updateData(linesDict, hostUpdateFile)
    print(len(newLinesDict), len(newLinesDict['Homo_sapiens']))
    return newLinesDict, cleanLabels