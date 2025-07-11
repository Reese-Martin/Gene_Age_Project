"""
This script loads the data frame generated in 4_Generate_DataFrame.py and then plots the distribution of 
Gene ages as a histogram

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter


def getcounts(df):
    c = Counter(df['AgeBin'].to_list())
    tmp = np.zeros([len(c), 2])
    tmp[:, 0] = list(c.keys())
    tmp[:, 1] = list(c.values())
    return tmp


save_folder = 'Images Combined/'

Cdata_folder = 'Data Celeg/'
Hdata_folder = 'Data Human/'
Mdata_folder = 'Data Mouse/'
Ddata_folder = 'Data Dmel/'
Adata_folder = 'Data Athal/'
Drdata_folder = 'Data Dreri/'
Ydata_folder = 'Data Yeast/'
image_folder = 'Images Combined/'
# load Gene age DF
Mdf = pd.read_pickle(Mdata_folder+"Gene_Age_DataFrame.pkl")
Hdf = pd.read_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")
Cdf = pd.read_pickle(Cdata_folder+"Gene_Age_DataFrame.pkl")
Ddf = pd.read_pickle(Ddata_folder+"Gene_Age_DataFrame.pkl")
Adf = pd.read_pickle(Adata_folder+"Gene_Age_DataFrame.pkl")
Drdf = pd.read_pickle(Drdata_folder+"Gene_Age_DataFrame.pkl")
Ydf = pd.read_pickle(Ydata_folder+"Gene_Age_DataFrame.pkl")
# load Gene age DF
Hdata = getcounts(Hdf)
Mdata = getcounts(Mdf)
Ddata = getcounts(Ddf)
Cdata = getcounts(Cdf)
Adata = getcounts(Adf)
Drdata = getcounts(Drdf)
Ydata = getcounts(Ydf)

Data = [Hdata, Mdata, Ddata, Cdata, Adata, Drdata]
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(12, 12))
cntr = 0
for ax in axes.flat:
    ax.bar(Data[cntr][:, 0], Data[cntr][:, 1])
    cntr = cntr + 1
# fig.tight_layout()
# set labels
# plt.subplots_adjust(left=0.1, bottom=0.2, right=.8, top=0.9, wspace=0, hspace=0)
axes[0, 0].set_title('Human')
axes[0, 1].set_title('Mouse')
axes[1, 0].set_title('Dmel')
axes[1, 1].set_title('C. elegans')
axes[2, 0].set_title('A. thaliana')
axes[2, 1].set_title('D. rerio')
plt.setp(axes[-1, :], xlabel='Age Bin')
plt.setp(axes[0, 0], ylabel='Biological Process Count')
plt.setp(axes[1, 0], ylabel='Biological Process Count')
plt.setp(axes[2, 0], ylabel='Biological Process Count')
fig.show()
fig.savefig(save_folder+"Gene_Age_Group_vs_Count.png")

