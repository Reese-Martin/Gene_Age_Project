#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: martinra4

this script uses the unique species from 2_Collect_Unique_Species.py (after processing
by TNRS: https://tree.opentreeoflife.org/curator/tnrs/) to generate a synthetic tree from the open tree of life and then retrieves the ages of
genes

"""
# trim the open tree of life super tree to only include the unq species we have
from opentree import OT
data_folder = 'Data Yeast/'

Rectified_Species = open(data_folder+'main.csv').readlines()
ott_ids = set()

for lin in Rectified_Species[1:]:  # skip the header
    lii = lin.split(',')  # split on commas
    try:
        ott_id = int(lii[2])  # grab the opentree id
        ott_ids.add(ott_id)  # add to the set
    except:
        pass

# use this try to capture any unknown query ids which can then be removed from the species list
try:
    treefile = data_folder+"Unq_Species_Tree.tre"
    # Get the synthetic tree from OpenTree
    output = OT.synth_induced_tree(ott_ids=list(ott_ids), label_format='name')
    output.tree.write(path=treefile, schema="newick")
    output.tree.print_plot(width=100)

except Exception as e:  # work on python 3.x
    Args = e.args
    spec_to_rem = list(Args)[0].split('\n')[2:-1]  # start at 2 to skip the opening lines of the eror message
    spec_to_rem = [s.strip('ott') for s in spec_to_rem]
    spec_to_rem = [s.strip(' ott') for s in spec_to_rem]
    for i in spec_to_rem:
        if int(i) in ott_ids:
            ott_ids.remove(int(i))
            # after problematic otts are trimmed retry the tree induction
    treefile = data_folder+"Unq_Species_Tree.tre"
    # Get the synthetic tree from OpenTree
    output = OT.synth_induced_tree(ott_ids=list(ott_ids), label_format='name')
    output.tree.write(path=treefile, schema="newick")
    output.tree.print_plot(width=100)

