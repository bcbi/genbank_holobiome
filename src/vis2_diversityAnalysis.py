"""
Diversity Analysis using Hill Numbers from Chao, et al. 2016 
for incidence data
Hill Numbers: q=0, q=1, and q=2 used 

Inputs: 
    1) PAUP Filename (ex: gi_phyla.nex)
    2) Host Updates CSV (ex. hostUpdates.csv)
    3) Output File name (full pathway with .png)

Author: Vivek Ramanan
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import seaborn as sns
import vis0_paupProcessing as paupProcessing

def qCalc(data, labels):
    """
    Calculates hill numbers for a given sample 
    @param data: dictionary from runningHillNumbers analysisDict, keys:hosts, values:lists of 0/1 data
    @param labels: column labels for given data
    @return q0,q1,q2 values all floats
    """
    q0=0
    q1=1
    q2=2
    total = np.zeros(len(labels))
    for host in data: 
        temp = np.array(data[host])
        total = np.add(temp, total)

    pi = []
    for k in range(len(total)):
        if total[k] != 0:
            pi.append(total[k] / len(data))

    richness= 0
    relativepi = []
    for i in range(len(pi)):
        calc = pi[i] / sum(pi)
        relativepi.append(calc)
        richness += calc**q0
    
    shannonAcc = 0
    simpsonAcc = 0
    for j in range(len(relativepi)):
        shanTemp = relativepi[j] * np.log(relativepi[j])
        shannonAcc += shanTemp
        simpTemp = relativepi[j]**2
        simpsonAcc += simpTemp
    shannon = np.exp(-1 * shannonAcc)
    simpson = 1 / (simpsonAcc)
    
    return richness, shannon, simpson

def seabornFigure(x,q0,q1,q2,outputFile):
    """
    Plots the figure using seaborn package
    @param x,q0,q1,q2: returned lists from runningHillNumbers function
    @param outputFile: FULL pathway to output file
    """
    sns.set(rc={'figure.figsize':(7,5)}, font_scale=3)
    sns.set_style("dark")
    sns.set_context("paper")
    sns.lineplot(x=xP, y=q0P, label="q=0 Richness")
    sns.lineplot(x=xP, y=q1P, label="q=1 Shannon")
    sns.lineplot(x=xP, y=q2P, label="q=2 Simpson")
    plt.savefig(outputFile, bbox_inches="tight")

def runningHillNumbers(data, labels):
    """
    Runs hill number process
    @param data: dictionary of PAUP data, keys are hosts and values are arrays of 0/1 data
    @param labels: list of column labels from PAUP data
    @return x,q0,q1,q2: all lists with float values
    """
    x = []
    q0 = []
    q1 = []
    q2 = []
    order = [] # if interested in order of which hosts were added to sample, print this
    analysisDict = {}
    for i in range(1, len(data)+1):
        key = random.choice(list(data.keys()))
        order.append(key)
        analysisDict[key] = data.pop(key)
        rich, shannon, simpson = qCalc(analysisDict, labels)
        x.append(i)
        q0.append(rich)
        q1.append(shannon)
        q2.append(simpson)
    return x,q0,q1,q2

def main():
    paupFile = sys.argv[1]
    hostUpdateFile = sys.argv[2]
    outputFile = sys.argv[3]

    data, labels = paupProcessing.PAUPprocess(paupFile, hostUpdateFile)
    x,q0,q1,q2 = runningHillNumbers(data, labels)
    seabornFigure(x,q0,q1,q2,outputFile)

main()