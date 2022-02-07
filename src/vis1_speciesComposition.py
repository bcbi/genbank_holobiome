"""
Visualization of Species Composition of GI Data

Inputs:
    1) PAUP Filename (ex: gi_species.nex)
    2) Host Updates CSV (ex. hostUpdates.csv)
    3) Taxonomy Dictionary text file (ex. taxonomyDict.txt)
    4) Output Filename (ex. should be png)

Author: Vivek Ramanan
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import vis0_paupProcessing as paupProcessing
from matplotlib.colors import ListedColormap

def createTaxDict(taxFile):
    """
    Based on taxonomy dictionary file created from speciesClean scripts in Cleaning Scripts
    Use given taxonomyDict.txt file unless recreating
    @return taxonomyDictionary as a dict
    """
    taxonomyDict = {}
    with open(taxFile, 'r') as f:
        for line in f:
            split = line.rstrip().split("\t")
            key = split[0]
            vals = split[1:]
            if "Sendai virus" in key:
                #print(key,flush=True)
                vals = ["genus:Respirovirus","family:Paramyxoviridae","order:Mononegavirales",\
                        "class:Monjiviricetes","phylum:Negarnaviricota","resolution:genus"]
            elif "Bastrovirus" in key:
                #print(key,flush=True)
                vals = ["family:Astroviridae","order:Stellavirales","class:Stelpaviricetes",\
                        "phylum:Pisuviricota","resolution:family"]
            keyDict = {}
            for value in vals:
                splitVal = value.split(":")
                keyDict[splitVal[0]] = splitVal[1]
            taxonomyDict[key] = keyDict
    print(len(taxonomyDict))
    return taxonomyDict

def valueCountsSpecies(newDF, l, tax):
    """
    Puts together the value counts for the given dataframe based on the hosts of interest
    and the taxonomyDict to get the phyla level grouping as well
    Percentages under 2% are grouped into the Other category
    @param newDF: species level dataframe
    @param l: group of interest of hosts
    @param tax: taxonomyDict
    @return vals (dict), otherMakeup(whats in the Other grouping), totalVals (final values)  
    """
    df = newDF[newDF.index.isin(l)]
    df.loc["Total"] = df.sum()
    vals = {}
    totalNum = 0
    for i in range(len(df.loc["Total"])):
        if df.loc["Total"][i] > 0:
            spec = df.columns[i]
            spec = spec.replace("_", " ")
            if spec in tax: 
                try: 
                    phylum = tax[spec]["phylum"]
                except:
                    continue
                    #print(spec, tax[spec])
                #print(spec, phylum)
                if phylum in vals: 
                    vals[phylum] += 1
                else:
                    vals[phylum] = 1
                totalNum += 1
    totalVals = {}
    otherVal = 0
    otherMakeup = []
    for key in vals.copy():
        percentage = 100*(float(vals[key]) / float(totalNum))
        if percentage <= 2.0:
            otherVal += vals[key]
            otherMakeup.append(key)
            vals.pop(key, None)
        totalVals[key] = percentage
    vals["Other"] = otherVal
    
    return vals, otherMakeup, totalVals

def createChart(cladeGroup, data, taxonomyDict, outputFile):
    """
    Creates final stacked 100 percent bar chart based on Seaborn packages and Pandas
    @param cladeGroup: clade dictionary keys:host group (str), values: arrays of host names
    @param data: species level dataframe
    @param taxonomyDict
    @param outputFile: FULL pathway to save chart at
    """
    dfData = []
    for clade in cladeGroup: 
        temp, other, totalTemp = valueCountsSpecies(data, cladeGroup[clade], taxonomyDict)
        relativeTemp = {}
        for val in temp:
            relativeTemp[val] = (temp[val] / sum(list(temp.values())))*100
        dfData.append(relativeTemp)

    tempDF = pd.DataFrame(dfData, index=list(cladeGroup.keys()))
    tempDF = tempDF.fillna(0)

    # Plotting
    sns.set(rc={'figure.figsize':(20,15)}, font_scale=2)
    ax = tempDF.plot(kind="bar", stacked=True, colormap=ListedColormap(sns.color_palette("twilight", 12)), rot=0)
    for rect in ax.patches:
        # Find where everything is located
        height = rect.get_height()
        width = rect.get_width()
        x = rect.get_x()
        y = rect.get_y()
        
        # The height of the bar is the data value and can be used as the label
        label_text = f'{height:.2f}%'  # f'{width:.2f}' to format decimal values
        
        # ax.text(x, y, text)
        label_x = x + width / 2
        label_y = y + height / 2
        
        # only plot labels greater than given width
        if height > 0.00:
            ax.text(label_x, label_y, label_text, ha='center', va='center', fontsize=20, color="w")

    plt.legend(loc="center right", bbox_to_anchor=(1.25, 0.5), ncol=1)
    plt.savefig(outputFile, bbox_inches="tight")
    plt.show()
    return

def processHostData(hostUpdateFile, data, taxonomyDict, outputFile):
    """
    Process the hostUpdateFile into groupings of interest for chart creation
    @param hostUpdateFile: FULL pathway to host Update file (same used on vis2)
    @param data: dataframe from paupProcess
    @param taxonomyDict
    @param outputFile: FULL pathway to saving chart
    """
    hostDF = pd.read_csv(hostUpdateFile)

    mammals = list(hostDF[hostDF['class'] == "Mammalia"]['Newick_label'])
    fish = list(hostDF[hostDF['class'] == "Actinopteri"]['Newick_label'])
    birds = list(hostDF[hostDF['class'] == "Aves"]['Newick_label'])
    bats = list(hostDF[hostDF['order'] == "Chiroptera"]['Newick_label'])
    giTotal = list(data.index.values[:-1])

    # Additional groupings I made based on interest
    primates = ["Pan_troglodytes", "Gorilla_gorilla", "Pan_paniscus", "Macaca_syvanus", "Papio_hamadryas", "Mandrillus_sphinx"]
    bats = ["Myotis_nattereri", "Myotis_daubentonii", "Plecotus_auritus", "Nyctalus_noctula", "Miniopterus_schreibersii", \
        "Pipistrellus_pygmaeus", "Eptesicus_serotinus", "Hipposideros_armiger", "Taphozous_perforatus"]
    domesticated = ["Homo_sapiens", "Sus_scrofa", "Bos_taurus", "Mus_musculus", "Ovis_aries", "Capra_aegrarus_hircus", \
               "Felis_catus", "Equus_caballus", "Canis_familiaris", "Gallus_gallus"]
    humans = ["Homo_sapiens"]
    
    # customize this based on choices 
    gi_clades = {"Overall":giTotal, "Mammals":mammals,  "Primates":primates, "Humans": humans, "Bats":bats}
    createChart(gi_clades, data, taxonomyDict, outputFile)
    return

def main():
    paupFile = sys.argv[1] # species level
    hostUpdateFile = sys.argv[2]
    data, labels = paupProcessing.PAUPprocess(paupFile, hostUpdateFile)
    taxonomyDict = createTaxDict(sys.argv[3])
    processHostData(hostUpdateFile, data, taxonomyDict, sys.argv[4])

main()
