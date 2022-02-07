"""
Network clustering using Markov Clustering and analysis of cluster results as heatmap
Additional presence-absence networks and mammalian networks produced

Inputs: 
1) 
    1) PAUP Filename (ex: gi_phyla.nex)
    2) Host Updates CSV (ex. hostUpdates.csv)
    3) Output DIRECTORY Name (ex. /users/name/Figures/ )

Author: Vivek Ramanan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import paupProcessing
from sklearn.metrics.pairwise import cosine_similarity
import itertools
import markov_clustering as mc
import networkx as nx

phylaComp = {'Bacteria': ['Synergistetes','Verrucomicrobia','Elusimicrobia','Tenericutes','Fusobacteria',
                'Cyanobacteria','Firmicutes','Lentisphaerae','Euryarchaeota','Spirochaetes','Fibrobacteres','Proteobacteria',
                'Chlamydiae','Bacteroidetes','Actinobacteria','Acidobacteria','Planctomycetes','Deferribacteres', 'DeinococcusThermus'],
                'Protista': ['Evosea','Cercozoa','Heterolobosea','Apicomplexa','Discosea',
                'Tubulinea','Parabasalia','Fornicata','Ciliophora','Euglenozoa'],
                'Invertebrates': ['Nematoda', 'Annelida', 'Acanthocephala'],
                'Viruses': ['Statovirus','Peploviricota','Preplasmiviricota','Iotatorquevirus','Anellovirus','Aquabirnavirus',
                'Negarnaviricota','Duplornaviricota','Artverviricota','Gyrovirus','Cressdnaviricota','Anelloviridae',
                'Cossaviricota','Kitrinoviricota','Uroviricota','Thetatorquevirus','Pisuviricota'],
                'Plants and Fungi': ['Blastocystis','Ascomycota','Dinophyceae','Chytridiomycota','Mucoromycota',
                'Basidiomycota','Chlorophyta','Streptophyta','Microsporidia'],
                'Phages': ['Phixviricota']}

def markovClustering(cosValues, hostDF, outputDir):
    """
    Runs Markov Clustering based on markov_clustering package
    @param cosValues: tuple dictionary with cosine similarity value for each pair
    @param hostDF: host dataframe
    @param outputDir: FULL pathway directory for creating new files
    @return G (network object), clusters (list of clusters)
    """
    G = nx.Graph()
    for b in cosValues: 
        if cosValues[b][0] > 0.66:
            G.add_edge(*b)
    matrix = nx.to_scipy_sparse_matrix(G)
    result = mc.run_mcl(matrix)           # run MCL with default parameters
    clusters = mc.get_clusters(result)    # get clusters

    nodes = list(G.nodes)
    for i in range(len(nodes)):
        for j in range(len(clusters)):
            if i in clusters[j]:
                G.nodes[nodes[i]]["cluster"] = j
                G.nodes[nodes[i]]["HostClass"] = hostDF.loc[hostDF['updatedScientificName'] == nodes[i]]['class'].values[0]
    
    nx.write_gml(G, outputDir+"markovClusteringNetwork.gml")
    return G, clusters

def cosineSimilarity(data, hostDF):
    """
    Find cosine similarity for all pairs of hosts using itertools
    @param data: data dictionary from readNextFile
    @param hostDF: host dataframe
    @return cosValues: dictionary with tuple pairs as keys, cosine similarity float as key
            cosMammals: same as cosValues but only with mammal hosts
    """
    cosValues = {}
    for d in itertools.combinations(list(data.keys()),2):
        temp = cosine_similarity([data[d[0]]], [data[d[1]]])[0]
        cosValues[d] = temp
        #if "Homo_sapiens" in d and temp[0] > 0.02:
            #print(d, temp)
    
    betaMammals = {}
    for b in data:
        if hostDF.loc[hostDF['updatedScientificName'] == b]['class'].values[0] == "Mammalia":
            betaMammals[b] = data[b]

    cosMammals = {}
    for d in itertools.combinations(list(betaMammals.keys()),2):
        tempBeta = cosine_similarity([betaMammals[d[0]]], [betaMammals[d[1]]])[0]
        cosMammals[d] = tempBeta
    return cosValues, cosMammals

def analyzeMakeup(nodes, data, labels):
    """
    Analyzes percentage makeup of the MCL clusters in groups
    @param nodes: nodes of the cluster
    @param data: data dict
    @param labels: all labels (list)
    @return totals: array of percentages per label (phyla)
    """
    totals = []
    for i in range(len(labels)):
        labelTotal = []
        for j in range(len(nodes)):
            if data[nodes[j]][i] == 1:
                labelTotal.append(1)
            else:
                labelTotal.append(0)
        percentage = np.mean(labelTotal)
        if percentage > 0.5: 
            print(labels[i], percentage)
        totals.append(percentage)
    return totals

def MCLAnalysis(G, clusters, hostDF, data, labels):
    """
    Analyzes MCL results across clusters to find similarities in microbiome composition
    @param G: network object of networkx
    @param clusters: lists of clusters
    @param hostDF: host dataframe
    @param data: data dict
    @param labels: labels list
    @return clusterPercentages: dictionary with keys as cluster number (starting from 1)
            and values being the percentage array across all phyla labels
    """
    nodes = list(G.nodes)
    clusterPercentages = {}
    for j in range(len(clusters)):
        print("Cluster",j)
        clusterNodes = []
        clusterClass = []
        clusterOrder = []
        for val in clusters[j]:
            tempNode = nodes[val]
            tempClass = hostDF.loc[hostDF['updatedScientificName'] == tempNode]['class'].values[0]
            tempOrder = hostDF.loc[hostDF['updatedScientificName'] == tempNode]['order'].values[0]
            clusterNodes.append(tempNode)
            clusterClass.append(tempClass)
            clusterOrder.append(tempOrder)
        if len(set(clusterClass)) == 1: 
            print("All %s" %(clusterClass[0]))
            print(pd.value_counts(pd.Series(clusterOrder)))
        else: 
            print(pd.value_counts(pd.Series(clusterClass)))
        percentages = analyzeMakeup(clusterNodes, data, labels)
        clusterPercentages[j+1] = percentages
        print()
    return clusterPercentages

def heatMap(clusterPercentages, labels, outputDir):
    """
    Creates heatmap based on MCLAnalysis output
    @param clusterPercentages: output of MCLAnalysis
    @param labels: labels list
    @param outputDir: FULL output directory pathway 
    """
    orderedColumns = []
    for val in list(phylaComp.values()):
        orderedColumns += val
    clusterDF = pd.DataFrame(clusterPercentages).T
    clusterDF.columns = labels
    clusterDF = clusterDF[orderedColumns]
    
    # Plotting HeatMap
    sns.set(rc={'figure.figsize':(30,10)}, font_scale=2)
    ax = sns.heatmap(clusterDF, cmap="viridis")
    ax.set( xlabel = "Microorganism Phyla", ylabel = "Cluster Number")
    fig = ax.get_figure()
    fig.savefig(outputDir+"heatmap.png", bbox_inches='tight') 
    return

def presenceAbsenceNetwork(data, labels, hostDF, outputDir):
    """
    Creates a presence absence network based on entire dataset, with taxonomy of hosts
    and species included for Cytoscape mapping and grouping
    @param data: data dict
    @param labels: labels list
    @param hostDF: host dataframe
    @param outputDir: FULL output directory for figures 
    """
    G = nx.Graph()
    for host in data:
        if host in hostDF['Newick_label'].values:
            for i in range(len(data[host])):
                if data[host][i] == 1.0:
                    G.add_edge(*(host, labels[i]))

    nodes = list(G.nodes)
    for i in range(len(nodes)):
        if nodes[i] in labels: 
            for key in phylaComp:
                if nodes[i] in phylaComp[key]: break
            G.nodes[nodes[i]]["Taxonomy"] = key
        else: 
            G.nodes[nodes[i]]["Taxonomy"] = hostDF.loc[hostDF['Newick_label'] == nodes[i]]['class'].values[0]

    nx.write_gml(G, outputDir+'PANetwork.gml')

def mammalNetwork(cosMammals,hostDF, outputDir):
    """
    Creates cosine similarity based network for mammalian hosts
    @param cosMammals: dictionary with tuple pairs of all hosts (mammals here) and cosine similarity as value
    @param hostDF: host dataframe
    @param outputDir: FULL output directory for figures 
    """
    G = nx.Graph()
    for b in cosMammals: 
        if cosMammals[b] > 0.66:
            G.add_edge(*b, weight=cosMammals[b][0])     
    nodes = list(G.nodes)
    for i in range(len(nodes)):
        G.nodes[nodes[i]]["HostOrder"] = hostDF.loc[hostDF['updatedScientificName'] == nodes[i]]['order'].values[0]     
    nx.write_gml(G, outputDir+'mammalNetwork.gml')
    return

def readNexFile(filename):
    """
    Read in process for Nexus file here - done on phylum level
    Slightly different from vis1/2 because does not turn into pandas dataframe for ease of analysis as a dict
    @param filename: full pathway to filename 
    @return data: dictionary of data with keys as host name, values as array of 0/1 data
            labels: list of column labels
    """
    labels = None
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            if "charlabels" in line:
                labels = line.replace(";","").rstrip().split(" ")[1:]
            elif "dimensions" in line:
                print(line)
                continue
            else:
                if "0" in line or "1" in line:
                    splitLine = line.rstrip().split(" ")
    return data, labels

def main():
    outputDir = sys.argv[3]
    hostDF = pd.read_csv(sys.argv[2])
    data, labels = readNexFile(sys.argv[1])
    cosValues, cosMammals = cosineSimilarity(data, hostDF)

    # General Networks
    presenceAbsenceNetwork(data, labels, hostDF, outputDir)
    mammalNetwork(cosMammals,hostDF, outputDir)

    # Clustered Network and Heatmap
    G, clusters = markovClustering(cosValues, hostDF, outputDir)
    clusterPercentages = MCLAnalysis(G, clusters, hostDF, data, labels)
    heatMap(clusterPercentages, labels, outputDir)

main()