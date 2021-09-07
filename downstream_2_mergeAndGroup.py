"""
Downstream 2 - Merges the isolation and tissue CSVs created by Downstream 1
and also adds the grouping category to the final CSV

Inputs: 1) isolation_BacteriaPreferred.csv
        2) tissue_bacteriaPreferred.csv
        3) sourceGroupings.txt
"""

import pandas as pd
import sys

def updateGroupings(df, groupingFile): 
    """
    @param df: the dataframe to be updated
    @param groupingFile: the file that has the organ groupings 
    @return df
    """

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
    """
    Combines the isolation and tissue CSVs
    @param file1: usually isolation
    @param file2: tissue
    @param groupingFile: to update the groups
    """
    df1 = pd.read_csv(file1, low_memory=False)
    df2 = pd.read_csv(file2, low_memory=False)
    df1['DataType'] = 'isolation'
    df2['DataType'] = 'tissue'

    merged = pd.concat([df1, df2])
    merged_Groups = updateGroupings(merged, groupingFile)

    return merged_Groups

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    groupingFile = sys.argv[3]
    merged = combineDatabases(file1, file2, groupingFile)
    merged.to_csv("merged_Host_Bacteria.csv",index=False)


main()
