#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: martinra4

this script takes in the trimmed oma-groups and generates a unique list of all
species present
"""

import csv
from omadb import Client
common_folder = 'Common Data/'
data_folder = 'Data Yeast/'

c = Client()  # create the oma client that will retrieve the species name for entered oma IDs

with open(data_folder+"YEAST_Containing_Groups.txt") as f:  # uses the trimmed oma-groups file
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)
# uniprot species stores the unique 5 letter entries for each omaID from D
UniProt_Species = []
Species_Protein = []
for Group in d[3:]:
    for entry in Group[2:]:
        Spec = entry[0:5]
        if Spec not in UniProt_Species:
            UniProt_Species.append(Spec)

with open(common_folder+"Oma_Species.txt") as f:  # Oma_species is a text file generated from the OmaDB website
    reader = csv.reader(f, delimiter="\t")
    a = list(reader)

Species_Name = []
for ID in UniProt_Species:
    for spec in a:
        if ID == spec[0]:
            if len(spec[2].split()) == 2 and spec[2] not in Species_Name:
                Species_Name.append(spec[2])
            elif spec[2].split()[0] + ' ' + spec[2].split()[1] not in Species_Name:
                Species_Name.append(spec[2].split()[0] + ' ' + spec[2].split()[1])

# save list of unique species from the list of species names
with open(data_folder+'Unique_species_names.txt', 'w') as file:
    file.writelines('\t'.join(Species_Name) + '\n')

with open(data_folder+'Unique_uniprot_ids.txt', 'w') as file:
    file.writelines('\t'.join(UniProt_Species) + '\n')
