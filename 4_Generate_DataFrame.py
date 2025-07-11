#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: martinra4

this script collects all relevant data in 1-3 and puts it into a single data frame
"""

import pandas as pd
import csv
import dendropy
from collections import Counter
import numpy as np

common_folder = 'Common Data/'
data_folder = 'Data Human/'

node_filter_fn = lambda nd: nd in descendants

# read csv file
with open(data_folder+"Human_Containing_Groups.txt") as f:  # uses the trimmed oma-groups file
    reader = csv.reader(f, delimiter="\t")
    OmaData = list(reader)
with open(common_folder+"Oma_Species.txt") as f:  # uses the trimmed oma-groups file
    reader = csv.reader(f, delimiter="\t")
    OmaSpecies = list(reader)

UniSpecies = []
SpeciesName = []
for i in OmaSpecies[3:]:
    UniSpecies.append(i[0])
    if len(i[2].split()) == 2:
        SpeciesName.append(i[2])
    elif i[2].split()[0] + ' ' + i[2].split()[1]:
        SpeciesName.append(i[2].split()[0] + ' ' + i[2].split()[1])

GroupNums = []
Fingerprints = []
Entry_IDs = []
for i in OmaData[3:]:
    GroupNums.append(i[0])
    Fingerprints.append(i[1])
    Entry_IDs.append(i[2:])

data = {'Group_ID': GroupNums, 'FingerPrints': Fingerprints, 'Entry_IDs': Entry_IDs}
df = pd.DataFrame(data=data)

EntrySpecies = []
for i in range(0, len(df)):
    if i % 100 == 0:
        print(i)
    tmp = df['Entry_IDs'][i]
    tmp2 = []
    for j in tmp:
        tmp2.append(SpeciesName[UniSpecies.index(j[0:5])])
    EntrySpecies.append(tmp2)

df['Species_Name'] = EntrySpecies

# load tree generated in 3_Generate_Tree.py
tree = dendropy.Tree.get(path=data_folder+"Unq_Species_Tree.tre", schema="newick")

# relational tree has no edge length so
# each edge between nodes gets a length of 1 for calculation of Node Distance from root
for edge in tree.postorder_edge_iter():
    edge.length = 1.0
# calculate distance from the root node for all nodes in the tree
tree.calc_node_root_distances()

# find the distance from the root node for each species a given gene is present in
All_Dist = []
for j in range(0, len(df['Species_Name'])):
    if j % 100 == 0:
        print(j)
    Dist = []
    for i in df['Species_Name'][j]:
        try:
            a = tree.mrca(taxon_labels=["Arabidopsis thaliana", i])
            # MRCAs.append(a)
            Dist.append(a.distance_from_root())
        except:
            Dist.append(100.0)
    All_Dist.append(Dist)

# find the minimum distance from the root for each ortholog
Min_Dist = []
for i in range(0, len(All_Dist)):
    Min_Dist.append(min(All_Dist[i]))
df['Root_Distances'] = All_Dist
df['Min_Dist'] = Min_Dist

# this section assigns an age bin to each protein so that results can be binned to avoid sparse samples skewing results

Age_Groups = df['Min_Dist'].to_list()
c = Counter(Age_Groups)
ages = list(c.keys())
ages.sort()
AgeBins = pd.DataFrame(columns=["Ages", "Counts"])
AgeBins["Ages"] = ages
AgeBins["Counts"] = [c[i] for i in ages]

bins = np.zeros(len(ages))
tmp = 0
bin = 0
for i in range(0, len(ages)):
    tmp = tmp + c[ages[i]]
    if tmp < 1000:
        bins[i] = bin
    if tmp > 1000:
        bins[i] = bin
        bin = bin +1
        tmp = 0

AgeBins['Bins'] = bins
df["AgeBin"] = [AgeBins[AgeBins['Ages'] == i]["Bins"].values[0] for i in df['Min_Dist']]

bins = np.zeros(len(ages))
tmp = 0
bin = 0
for i in range(0, len(ages)):
    tmp = tmp + c[ages[i]]
    if tmp < 500:
        bins[i] = bin
    if tmp > 500:
        bins[i] = bin
        bin = bin +1
        tmp = 0

AgeBins['Bins'] = bins
df["AgeBin_500"] = [AgeBins[AgeBins['Ages'] == i]["Bins"].values[0] for i in df['Min_Dist']]
# for both Athaliana and Dmel the last set of genes is very small (<100 genes) so you need to manually
# combin that age bin with the second to last age bin.


# to trim to Humans use 'HUMAN'
# to trim to Mouse use 'MOUSE'
# to trim to Drosophila Melanogaster use "DROME"
# to trim to C. elegans use "CAEEL"
# to trim to C. elegans use "CAEEL"
# to trim to C. elegans use "CAEEL"
# to trim to A. thaliana use "ARATH"
# to trim to D. Rerio use DANRE
# to trim to S. cerivisiae use YEAST
GeneLabel = []
sub = 'ARATH'
for i in df['Entry_IDs']:
    GeneLabel.append(list(filter(lambda x: sub in x, i)))

df['Gene'] = GeneLabel
# Save DF as pickle
df.to_pickle(data_folder+"Gene_Age_DataFrame.pkl")
