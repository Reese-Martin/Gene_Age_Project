# this script loads the dataframes for each species, then for the common categories add columns to the dataframe
# indicating whether a gene participates in the given process
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


def Process_Counts(df, x):
    for i in x:
        tmp = [lst.count(i) for lst in df['Broad_BioProcesses']]
        df[i] = tmp
    return df


Cdata_folder = 'Data Celeg/'
Hdata_folder = 'Data Human/'
Mdata_folder = 'Data Mouse/'
Ddata_folder = 'Data Dmel/'
Adata_folder = 'Data Athal/'
Drdata_folder = 'Data Dreri/'
image_folder = 'Images Combined/'

# load Gene age DF
Mdf = pd.read_pickle(Mdata_folder+"Gene_Age_DataFrame.pkl")
Hdf = pd.read_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")
Cdf = pd.read_pickle(Cdata_folder+"Gene_Age_DataFrame.pkl")
Ddf = pd.read_pickle(Ddata_folder+"Gene_Age_DataFrame.pkl")
Adf = pd.read_pickle(Adata_folder+"Gene_Age_DataFrame.pkl")
Drdf = pd.read_pickle(Drdata_folder+"Gene_Age_DataFrame.pkl")

ABPs = list(set([i for sublist in Adf['Broad_BioProcesses'] for i in sublist]))
DBPs = list(set([i for sublist in Ddf['Broad_BioProcesses'] for i in sublist]))
DrBPs = list(set([i for sublist in Drdf['Broad_BioProcesses'] for i in sublist]))
CBPs = list(set([i for sublist in Cdf['Broad_BioProcesses'] for i in sublist]))
MBPs = list(set([i for sublist in Mdf['Broad_BioProcesses'] for i in sublist]))
HBPs = list(set([i for sublist in Hdf['Broad_BioProcesses'] for i in sublist]))

# get the common processes for each species
collection = [ABPs, DBPs, DrBPs, MBPs, CBPs, HBPs]
collection = [i for sublist in collection for i in sublist]
counts = Counter(collection)
Common = []
for i in counts.keys():
    if counts[i] == 6:  # change value to the number of data sources compared
        Common.append(i)

Adf = Process_Counts(Adf, Common)
Cdf = Process_Counts(Cdf, Common)
Ddf = Process_Counts(Ddf, Common)
Drdf = Process_Counts(Drdf, Common)
Hdf = Process_Counts(Hdf, Common)
Mdf = Process_Counts(Mdf, Common)

Adf.to_pickle(Adata_folder+"Gene_Age_DataFrame.pkl")
Cdf.to_pickle(Cdata_folder+"Gene_Age_DataFrame.pkl")
Ddf.to_pickle(Ddata_folder+"Gene_Age_DataFrame.pkl")
Drdf.to_pickle(Drdata_folder+"Gene_Age_DataFrame.pkl")
Hdf.to_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")
Mdf.to_pickle(Mdata_folder+"Gene_Age_DataFrame.pkl")

# generate barplot of counts of genes that have a given association
Adata = sum(Adf[Common].values)
Cdata = sum(Cdf[Common].values)
Ddata = sum(Ddf[Common].values)
Drdata = sum(Drdf[Common].values)
Hdata = sum(Hdf[Common].values)
Mdata = sum(Mdf[Common].values)

Data = [Hdata, Mdata, Ddata, Cdata, Adata, Drdata]

fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(12, 20), sharex=True)
cntr = 0
for ax in axes.flat:
    ax.bar(range(1, len(Common)+1), Data[cntr])
    ax.set_ylim(0, 12500)
    ax.plot([0, len(Common)], [3000, 3000], '--')
    ax.plot([0, len(Common)], [2000, 2000], '--', color='black')
    if cntr in range(0, 4): #remove x ticks from top plots
        ax.tick_params(bottom=False)
    else:
        ax.set_xticks(np.arange(len(Common))+1, labels=Common)
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    cntr = cntr + 1
# fig.tight_layout()
axes[0, 0].set_title('Human')
axes[0, 1].set_title('Mouse')
axes[1, 0].set_title('Dmel')
axes[1, 1].set_title('C. elegans')
axes[2, 0].set_title('A. thaliana')
axes[2, 1].set_title('D. rerio')
plt.setp(axes[0, 0], ylabel='Biological Process Count')
plt.setp(axes[1, 0], ylabel='Biological Process Count')
fig.show()
fig.savefig(image_folder+"Count_Genes_in_Processes.pdf")
