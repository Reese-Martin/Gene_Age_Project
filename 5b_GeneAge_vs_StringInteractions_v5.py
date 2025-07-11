"""
This script loads the data frame generated in 4_Generate_DataFrame.py and then plots the distribution of
Gene ages against the number of biological processes
v4 cuts out unused code, and switches to log scale violin plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def AgeSortedArray(df):
    Age_PPI = np.zeros([len(df), 2])
    Age_PPI[:, 0] = np.asarray(list(df['AgeBin']))
    Age_PPI[:, 1] = np.asarray(list(df['MedConf_PPI'])) + np.asarray(list(df['HighConf_PPI']))
    # Age_PPI[:, 1] = np.asarray(df['UHC_PPI'])
    Age_PPI_sorted = Age_PPI[np.argsort(Age_PPI[:, 0])]
    Age_PPI_sorted[np.argwhere(Age_PPI_sorted[:, 1] == 0), 1] = np.nan
    return Age_PPI_sorted


def anova(x):
    tmp = x[~np.isnan(x[:, 1]), :]
    df = pd.DataFrame(columns=['Age', 'PPI'], data=tmp)
    lm = ols('PPI~Age', data=df).fit()
    table = sm.stats.anova_lm(lm)
    print(table)
    Age_pval = table['PR(>F)']['Age']
    if Age_pval < .05:
        a = pairwise_tukeyhsd(df['PPI'], df['Age']).summary().data
    else:
        a = 0
    return Age_pval, a


def logBinMedian(x, cols, bins):
    data_array = np.zeros([bins, cols])
    for i in range(0, bins):
        data_array[i, 0] = i
        for j in range(1, cols):
            data_array[i, j] = np.log(np.nanmedian(x[np.where(x[:, 0] == i)][:, j]))
    return data_array


def logviolins(x, bins):
    data_array = []
    for i in range(0, bins):
        tmp = x[np.where(x[:, 0] == i)][:, 1]
        tmp = tmp[~np.isnan(tmp)]
        data_array.append(list(np.log(tmp)))
    return data_array


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

# count the number of entries there are for each age group
HAge_Sorted = AgeSortedArray(Hdf)
Hpval, tmp = anova(HAge_Sorted)
HlogMedian_Age = logBinMedian(HAge_Sorted, 2, int(max(HAge_Sorted[:, 0])+1))
Hviolins = logviolins(HAge_Sorted, int(max(HAge_Sorted[:, 0])+1))

MAge_Sorted = AgeSortedArray(Mdf)
Mpval, tmp = anova(MAge_Sorted)
MlogMedian_Age = logBinMedian(MAge_Sorted, 2, int(max(MAge_Sorted[:, 0])+1))
Mviolins = logviolins(MAge_Sorted, int(max(MAge_Sorted[:, 0])+1))

DAge_Sorted = AgeSortedArray(Ddf)
Dpval, tmp = anova(DAge_Sorted)
DlogMedian_Age = logBinMedian(DAge_Sorted, 2, int(max(DAge_Sorted[:, 0])+1))
Dviolins = logviolins(DAge_Sorted, int(max(DAge_Sorted[:, 0])+1))

CAge_Sorted = AgeSortedArray(Cdf)
Cpval, tmp = anova(CAge_Sorted)
ClogMedian_Age = logBinMedian(CAge_Sorted, 2, int(max(CAge_Sorted[:, 0])+1))
Cviolins = logviolins(CAge_Sorted, int(max(CAge_Sorted[:, 0])+1))

AAge_Sorted = AgeSortedArray(Adf)
Apval, tmp = anova(AAge_Sorted)
AlogMedian_Age = logBinMedian(AAge_Sorted, 2, int(max(AAge_Sorted[:, 0])+1))
Aviolins = logviolins(AAge_Sorted, int(max(AAge_Sorted[:, 0])+1))

DrAge_Sorted = AgeSortedArray(Drdf)
Drpval, tmp = anova(DrAge_Sorted)
DrlogMedian_Age = logBinMedian(DrAge_Sorted, 2, int(max(DrAge_Sorted[:, 0])+1))
Drviolins = logviolins(DrAge_Sorted, int(max(DrAge_Sorted[:, 0])+1))

Data = [Hviolins, Mviolins, Drviolins, Dviolins, Cviolins, Aviolins]
pvals = [Hpval, Mpval, Drpval, Dpval, Cpval, Apval]
Meds = [HlogMedian_Age, MlogMedian_Age, DrlogMedian_Age, DlogMedian_Age, ClogMedian_Age, AlogMedian_Age]
xlims = [[0, 11], [0, 13], [0, 13], [0, 9], [0, 9], [0, 12.1]]

hticks = [1, 10]
hlabels = ['Oldest', 'Youngest']

mticks = [1, 12]
mlabels = ['Oldest', 'Youngest']

drticks = [1, 12]
drlabels = ['Oldest', 'Youngest']

dticks = [1, 8]
dlabels = ['Oldest', 'Youngest']

cticks = [1, 8]
clabels = ['Oldest', 'Youngest']

aticks = [1, 11]
alabels = ['Oldest', 'Youngest']

ticks = [hticks, mticks, drticks, dticks, cticks, aticks]
labels = [hlabels, mlabels, drlabels, dlabels, clabels, alabels]
# making violin plots of data
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(16, 8))
cntr = 0
for ax in axes.flat:
    parts = ax.violinplot(Data[cntr], showmeans=False, showmedians=False, showextrema=False)
    ax.plot(Meds[cntr][:, 0]+1, Meds[cntr][:, 1], c='b', linewidth=.5)
    # ax.set_ylim([0, 6.3])
    ax.set_xlim(xlims[cntr])
    ax.set_xticks(ticks[cntr], labels[cntr])
    cntr2 = 1
    for i in Data[cntr]:
        if i:
            quartile1, medians, quartile3 = np.percentile(i, [25, 50, 75])
            whiskers = np.array(adjacent_values(i, quartile1, quartile3))
            whiskersMin, whiskersMax = whiskers[0], whiskers[1]

            ax.scatter(cntr2, medians, marker='_', color='white', s=5, zorder=3)
            ax.vlines(cntr2, quartile1, quartile3, color='k', linestyle='-', lw=2.5)
            ax.vlines(cntr2, whiskersMin, whiskersMax, color='k', linestyle='-', lw=.5)
            cntr2 = cntr2+1
        else:
            print('skipped')
    cntr = cntr + 1

axes[0, 0].set_title('$H. sapiens$')
axes[0, 1].set_title('$M. musculus$')
axes[0, 2].set_title('$D. rerio$')
axes[1, 0].set_title('$D. melanogaster$')
axes[1, 1].set_title('$C. elegans$')
axes[1, 2].set_title('$A. thaliana$')
plt.setp(axes[-1, :], xlabel='Age Bin')
plt.setp(axes[0, 0], ylabel='log(PPI)')
plt.setp(axes[1, 0], ylabel='log(PPI)')
fig.show()
fig.savefig(image_folder+"S1_Age_vs_PPI_boxplot.tiff")
