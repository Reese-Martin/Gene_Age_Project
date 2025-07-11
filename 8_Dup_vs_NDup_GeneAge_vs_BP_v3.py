# this script re-runs the bioprocess analysis but separates the violins based on duplication status
"""
This script loads the data frame generated in 4_Generate_DataFrame.py and then plots the distribution of
Gene ages against the number of biological processes
v4 cuts out unused code, and switches to log scale violin plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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


def anova_2w(x, x2):
    tmp = np.vstack([x, x2])
    status = ['Dup']*len(x) + ['Single']*len(x2)
    data = [[tmp[i, 0], tmp[i, 1], status[i]] for i in range(len(tmp))]

    df = pd.DataFrame(columns=['Ages', 'BPs', 'status'], data=data)

    lm = ols('BPs~Ages+C(status)', data=df).fit()
    table = sm.stats.anova_lm(lm)
    print(table)
    Age_pval = table['PR(>F)'][['Ages', 'C(status)']].values

    return Age_pval


def D_AgeSortedArray(df, thresh):
    Age_BP = np.zeros([len(df[df['Paralogs'] > thresh]), 2])
    Age_BP[:, 0] = np.asarray([df.loc[i, 'AgeBin'] for i in df.index if df.loc[i, 'Paralogs'] > thresh])
    Age_BP[:, 1] = np.asarray([len(df.loc[i, 'Assc_BioProcesses']) for i in df.index if df.loc[i, 'Paralogs'] > thresh])
    Age_BP_sorted = Age_BP[np.argsort(Age_BP[:, 0])]
    Age_BP_sorted[np.argwhere(Age_BP_sorted[:, 1] == 0), 1] = np.nan
    return Age_BP_sorted


def logviolins(x, bins):
    data_array = []
    for i in range(0, bins):
        tmp = x[np.where(x[:, 0] == i)][:, 1]
        tmp = tmp[~np.isnan(tmp)]
        data_array.append(list(np.log(tmp)))
    return data_array


def logBinMedian(x, cols, bins):
    data_array = np.zeros([bins, cols])
    for i in range(0, bins):
        data_array[i, 0] = i
        for j in range(1, cols):
            data_array[i, j] = np.log(np.nanmedian(x[np.where(x[:, 0] == i)][:, j]))
    return data_array


def ND_AgeSortedArray(df, thresh):
    Age_BP = np.zeros([len(df[df['Paralogs'] <= thresh]), 2])
    Age_BP[:, 0] = np.asarray([df.loc[i, 'AgeBin'] for i in df.index if df.loc[i, 'Paralogs'] <= thresh])
    Age_BP[:, 1] = np.asarray([len(df.loc[i, 'Assc_BioProcesses']) for i in df.index if df.loc[i, 'Paralogs'] <= thresh])
    Age_BP_sorted = Age_BP[np.argsort(Age_BP[:, 0])]
    Age_BP_sorted[np.argwhere(Age_BP_sorted[:, 1] == 0), 1] = np.nan
    return Age_BP_sorted


Hdata_folder = 'Data Human/'
Mdata_folder = 'Data Mouse/'
Ddata_folder = 'Data Dmel/'
Drdata_folder = 'Data Dreri/'
image_folder = 'Images Combined/'
# load Gene age DF
Hdf = pd.read_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")
Mdf = pd.read_pickle(Mdata_folder+"Gene_Age_DataFrame.pkl")
Drdf = pd.read_pickle(Drdata_folder+"Gene_Age_DataFrame.pkl")
Ddf = pd.read_pickle(Ddata_folder+"Gene_Age_DataFrame.pkl")

# count the number of entries there are for each age group
thresh = 0
D_HAge_Sorted = D_AgeSortedArray(Hdf, thresh)
D_HlogMedian_Age = logBinMedian(D_HAge_Sorted, 2, int(max(D_HAge_Sorted[:, 0])+1))
D_Hviolins = logviolins(D_HAge_Sorted, int(max(D_HAge_Sorted[:, 0])+1))
ND_HAge_Sorted = ND_AgeSortedArray(Hdf, thresh)
ND_HlogMedian_Age = logBinMedian(ND_HAge_Sorted, 2, int(max(ND_HAge_Sorted[:, 0])+1))
ND_Hviolins = logviolins(ND_HAge_Sorted, int(max(ND_HAge_Sorted[:, 0])+1))
Hpvals = anova_2w(D_HAge_Sorted, ND_HAge_Sorted)

D_Mage_Sorted = D_AgeSortedArray(Mdf, thresh)
D_MlogMedian_Age = logBinMedian(D_Mage_Sorted, 2, int(max(D_Mage_Sorted[:, 0])+1))
D_Mviolins = logviolins(D_Mage_Sorted, int(max(D_Mage_Sorted[:, 0])+1))
ND_MAge_Sorted = ND_AgeSortedArray(Mdf, thresh)
ND_MlogMedian_Age = logBinMedian(ND_MAge_Sorted, 2, int(max(ND_MAge_Sorted[:, 0])+1))
ND_Mviolins = logviolins(ND_MAge_Sorted, int(max(ND_MAge_Sorted[:, 0])+1))
Mpvals = anova_2w(D_Mage_Sorted, ND_MAge_Sorted)

D_Dage_Sorted = D_AgeSortedArray(Ddf, thresh)
D_DlogMedian_Age = logBinMedian(D_Dage_Sorted, 2, int(max(D_Dage_Sorted[:, 0])+1))
D_Dviolins = logviolins(D_Dage_Sorted, int(max(D_Dage_Sorted[:, 0])+1))
ND_DAge_Sorted = ND_AgeSortedArray(Ddf, thresh)
ND_DlogMedian_Age = logBinMedian(ND_DAge_Sorted, 2, int(max(ND_DAge_Sorted[:, 0])+1))
ND_Dviolins = logviolins(ND_DAge_Sorted, int(max(ND_DAge_Sorted[:, 0])+1))
Dpvals = anova_2w(D_Dage_Sorted, ND_DAge_Sorted)

D_Drage_Sorted = D_AgeSortedArray(Drdf, thresh)
D_DrlogMedian_Age = logBinMedian(D_Drage_Sorted, 2, int(max(D_Drage_Sorted[:, 0]) + 1))
D_Drviolins = logviolins(D_Drage_Sorted, int(max(D_Drage_Sorted[:, 0]) + 1))
ND_DRAge_Sorted = ND_AgeSortedArray(Drdf, thresh)
ND_DrlogMedian_Age = logBinMedian(ND_DRAge_Sorted, 2, int(max(ND_DRAge_Sorted[:, 0]) + 1))
ND_Drviolins = logviolins(ND_DRAge_Sorted, int(max(ND_DRAge_Sorted[:, 0]) + 1))
Drpvals = anova_2w(D_Drage_Sorted, ND_DRAge_Sorted)

D_Data = [D_Hviolins, D_Mviolins, D_Drviolins, D_Dviolins]
ND_Data = [ND_Hviolins, ND_Mviolins, ND_Drviolins, ND_Dviolins]
pvals = [Hpvals, Mpvals, Drpvals, Dpvals]
D_Meds = [D_HlogMedian_Age, D_MlogMedian_Age, D_DrlogMedian_Age, D_DlogMedian_Age]
ND_Meds = [ND_HlogMedian_Age, ND_MlogMedian_Age, ND_DrlogMedian_Age, ND_DlogMedian_Age]


hticks = [1.25, 10.25]
hlabels = ['Oldest', 'Youngest']

mticks = [1.25, 12.25]
mlabels = ['Oldest', 'Youngest']

drticks = [1.25, 12.25]
drlabels = ['Oldest', 'Youngest']

dticks = [1.25, 8.25]
dlabels = ['Oldest', 'Youngest']


ticks = [hticks, mticks, drticks, dticks]
labels = [hlabels, mlabels, drlabels, dlabels]
# making violin plots of data
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12, 8))
cntr = 0
for ax in axes.flat:
    if cntr < 5:
        positions = np.asarray(list(range(len(D_Data[cntr]))))+1
        ND_parts = ax.violinplot(ND_Data[cntr], positions=positions, showmeans=False, showmedians=False, showextrema=False)
        D_parts = ax.violinplot(D_Data[cntr], positions=positions+.5, showmeans=False, showmedians=False, showextrema=False)

        ax.plot(ND_Meds[cntr][:, 0]+1, ND_Meds[cntr][:, 1], c='b', linewidth=.5)
        ax.plot(D_Meds[cntr][:, 0]+1.5, D_Meds[cntr][:, 1], c='orange', linewidth=.5)
        ax.set_ylim([0, 6.2])
        ax.set_xticks(ticks[cntr], labels[cntr])
        cntr2 = 1
        for i in ND_Data[cntr]:
            quartile1, medians, quartile3 = np.percentile(i, [25, 50, 75])
            whiskers = np.array(adjacent_values(i, quartile1, quartile3))
            whiskersMin, whiskersMax = whiskers[0], whiskers[1]

            ax.scatter(cntr2, medians, marker='_', color='white', s=5, zorder=3)
            ax.vlines(cntr2, quartile1, quartile3, color='k', linestyle='-', lw=2.5)
            ax.vlines(cntr2, whiskersMin, whiskersMax, color='k', linestyle='-', lw=.5)
            cntr2 = cntr2 + 1

        cntr2 = 1
        for i in D_Data[cntr]:
            quartile1, medians, quartile3 = np.percentile(i, [25, 50, 75])
            whiskers = np.array(adjacent_values(i, quartile1, quartile3))
            whiskersMin, whiskersMax = whiskers[0], whiskers[1]

            ax.scatter(cntr2+.5, medians, marker='_', color='white', s=5, zorder=3)
            ax.vlines(cntr2+.5, quartile1, quartile3, color='grey', linestyle='-', lw=2.5)
            ax.vlines(cntr2+.5, whiskersMin, whiskersMax, color='grey', linestyle='-', lw=.5)
            cntr2 = cntr2 + 1
        cntr = cntr + 1

# fig.delaxes(axes[1, 2])
axes[0, 0].set_title('$H. sapiens$')
axes[0, 1].set_title('$M. musculus$')
axes[1, 0].set_title('$D. rerio$')
axes[1, 1].set_title('$D. melanogaster$')
plt.setp(axes[-1, :], xlabel='Age Bin')
plt.setp(axes[0, 0], ylabel='log(Biological Process Count)')
plt.setp(axes[1, 0], ylabel='log(Biological Process Count)')
single = mpatches.Patch(color='blue', alpha=.3, label='Singletons')
dup = mpatches.Patch(color='orange', alpha=.3, label='Duplicates')
fig.legend(handles=[single, dup], loc=[.795, .4])
fig.show()
fig.savefig(image_folder+"F3_Age_vs_BP_paralogs_boxplot.tiff", dpi=300)

