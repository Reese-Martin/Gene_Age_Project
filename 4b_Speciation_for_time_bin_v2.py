# this script is used to extract species specific to each age bin. Then used to generate a age bin to taxonomic
# rank conversion

import pandas as pd
from collections import Counter

# for each age bin, get the largest and smallest dist in that bin
data_folder = 'Data Athal/'
df = pd.read_pickle(data_folder+'Gene_Age_DataFrame.pkl')
dist_range = []
bins = list(Counter(df['AgeBin']).keys())
bins.sort()
for i in bins:
    a = df[df['AgeBin'] == i]['Min_Dist'].to_list()
    lb = min(a)
    ub = max(a)
    if ub == 100:
        tmp = list(set(a))
        tmp.remove(100)
        ub = max(tmp)
    dist_range.append([lb, ub])
bin_ranges = list(zip(bins, dist_range))
age_df = pd.DataFrame(bin_ranges, columns=['bins', 'range'])
# for each distance, retrieve a species so that age can be determined using the https://timetree.org/
tmp_lb_Species = []
tmp_ub_Species = []
for i in bin_ranges:
    a = df[df['Min_Dist'] == i[1][0]]['Root_Distances'].to_list()
    b = df[df['Min_Dist'] == i[1][0]]['Species_Name'].to_list()
    tmp_Species = []
    for j in range(0, len(a)):
        for k in range(0, len(a[j])):
            if a[j][k] == i[1][0]:
                tmp_Species.append(b[j][k])
    tmp = Counter(tmp_Species)
    tmp_lb_Species.append(max(tmp, key=tmp.get))

    c = df[df['Min_Dist'] == i[1][1]]['Root_Distances'].to_list()
    d = df[df['Min_Dist'] == i[1][1]]['Species_Name'].to_list()
    tmp_Species = []
    for j in range(0, len(c)):
        for k in range(0, len(c[j])):
            if c[j][k] == i[1][1]:
                tmp_Species.append(d[j][k])
    tmp = Counter(tmp_Species)
    tmp_ub_Species.append(max(tmp, key=tmp.get))

age_df['Species'] = list(zip(tmp_lb_Species, tmp_ub_Species))
# taxonomic ranks were manually determined based on the taxonomy of species in a bin
# human taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Phylum: Chordata', 'Phylum: Chordata',
#          'Phylum: Chordata', 'Phylum: Chordata', 'Class: Mammalia', 'Family: Hominidae']
# mouse taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Phylum: Chordata', 'Phylum: Chordata',
#          'Phylum: Chordata', 'Phylum: Chordata', 'Class: Mammalia', 'Class: Mammalia', 'Family: Muridae']
# danio rerio taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Kingdom: Animalia',
#          'Phylum: Chordata', 'Phylum: Chordata', 'Phylum: Chordata', 'Class: Actinopterygii', 'Class: Actinopterygii', 'Famliy: Cyprinidae' ]
# Drosophila melanogaster taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Phylum: Arthropoda', 'Order: Diptera',
#          'Genus: Drosophila', 'Species Group: melanogaster']
# Caenorhabditis elegans taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Kingdom: Animalia', 'Kingdom: Animalia', 'Phylum: Nematoda', 'Class: Chromadorea',
#          'Genus: Caenorhabditis', 'Genus: Caenorhabditis']
# Arabidopsis thaliana taxonomic ranks
# ranks = ['All Life', 'Domain: Eukaryotes', 'Domain: Eukaryotes', 'Kingdom: Plantae', 'Phylum: tracheophyta', 'Clade: Angiosperms',
#          'Clade: Angiosperms', 'Clade: Eudicots', 'Clade: Rosids', 'Family: Brassicaceae', 'Genus: Arabidopsis']

age_df['Rank'] = ranks
age_df.to_pickle(data_folder+'Bins_to_True_age.pkl')
