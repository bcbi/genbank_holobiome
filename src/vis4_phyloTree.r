---
title: "Tree Visualization Script"
author: "Vivek Ramanan"
date: "3/15/2021"
output: html_document
---

## This script creates two phylogenetic trees: one based on microbiome analysis, the output of the downstream steps
## and the other based on the SAME hosts from the microbiome analysis but based on evolutionary distances from 
## TimeTree.org. To see how to create that tree, please refer to the README

```{r}
library(ggtree)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(ape)
library(dendextend)
library(dplyr)
```
## Reads in host data file: slightly different organization than hostUpdates.csv
## Includes evolutionary tree hosts, based on different naming conventions from TimeTree.org
```
d <- read.csv("hostDataforR.csv")
```

## Reading in the microbiome based tree
```{r}
tree <- treeio::read.nexus("<microbiomeGItree>.tre")
```
## Tree with branch lengths included
```{r fig1, fig.height=20, fig.width=10}
ggtree(tree) %<+% d + 
  geom_tippoint(aes(shape=class, color=class), size=4) + geom_tiplab(offset=.005)
ggsave("gi_branchLengths.png")
```
## Tree in circular layout, including colored tips based on taxonomical class and labels for host names
```{r fig2, fig.height=15, fig.width=15}
ggtree(tree, branch.length = "none", layout="circular") %<+% d + 
  geom_tippoint(aes(color=class), size=4) + geom_tiplab(offset=1) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text=element_text(size=14)) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("gi_circular.png")
```
## Tree in circular layout, including colored tips based on taxonomical class and with NO LABELS
```{r fig3, fig.height=15, fig.width=15}
ggtree(tree, branch.length = "none", layout="circular") %<+% d + 
  geom_tippoint(aes(color=class), size=4) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text=element_text(size=14)) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("gi_circular_nolabels.png")
```

## Reading in evolutionary tree from TIMETREE
```{r}
hosttree <- treeio::read.newick("<evolutionaryTree>.tre")
```
## Evolutionary tree, circular layout, colored tips from taxonomy, labels for host names
```{r fig4, fig.height=15, fig.width=15}
ggtree(hosttree, layout="circular") %<+% d+
  geom_tippoint(aes(color=class), size=4) + geom_tiplab(offset=15) + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text=element_text(size=14)) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("evo_circular.png")
```
## Evolutionary tree, circular layout, colored tips from taxonomy, NO LABELS
```{r fig5, fig.height=15, fig.width=15}
ggtree(hosttree, layout="circular") %<+% d+
  geom_tippoint(aes(color=class), size=4) + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text=element_text(size=14)) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("evo_circular_nolabels.png")
```
