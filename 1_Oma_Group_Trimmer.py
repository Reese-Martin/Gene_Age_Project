#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: martinra4

the purpose of this script is to load the OMA-Groups.txt file and then trim that 
file to include only the groups that include a Human protein. 

Note: the first five letters of each entry in the oma groups file are the UniProtKB
species identifiers
"""
import csv

data_folder = 'Common Data/'
save_folder = 'Data Dreri/'
with open(data_folder+"oma-groups.txt") as f:  # oma-groups is downloaded from the OMA database website
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)

# create a list and then append the first 3 rows to list to maintain sorting
# info in trimmed list of groups
Trimmed = []
Trimmed.append(d[0])
Trimmed.append(d[1])
Trimmed.append(d[2])

# if group contains a human protein (denoted with 'HUMAN') then append that group
# to trimmed
# to trim to Humans use 'HUMAN'
# to trim to Mouse use 'MOUSE'
# to trim to Drosophila Melanogaster use "DROME"
# to trim to C. elegans use "CAEEL"
# to trim to Danio Rerio use DANRE
# to trim to Arabidopsis thaliana use ARATH
# to trim to saccharomyces use YEAST
for group in d:
    tmp = '\t'.join(group)
    if 'DANRE' in tmp:
        Trimmed.append(group)

# save trimmed group list to same directory original file is located in
with open(save_folder+'YEAST_Containing_Groups.txt', 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in Trimmed)
