"""
Species Cleaning 1
Create microSpeciesNames.txt file and splits it accordingly as well
Author: Vivek Ramanan
"""

import pandas as pd
import sys, os

def splitSpecies(dfFile, outputName): 
    """
    Splits the list of species from the dataframe CSV
    @param dfFile: dataframe file
    @param outputName: name of the output file
    """
    df = pd.read_csv(dfFile, low_memory=False, usecols=['species','host_cleaned'])
    df = df.dropna(subset=['host_cleaned'])
    speciesNames = list(df['species'].value_counts().index)
    with open(outputName, 'w') as f: 
        for name in speciesNames: 
            f.write(name + "\n")
    
    # split file into groups of 15,000
    command = "split -a 1 -l 15000 " + outputName + " microSpecies"
    os.system(command)

def main():
    splitSpecies(sys.argv[1], sys.argv[2])

main()