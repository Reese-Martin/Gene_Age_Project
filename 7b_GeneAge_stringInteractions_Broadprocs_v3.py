# this script generates plots of the biological processes of a gene separated by the broad processes
# that gene participates in

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import random
import statsmodels.api as sm
from statsmodels.formula.api import ols


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def AgeSortedArray(df, trait):
    Age_PPI = np.zeros([len(df), 2])
    Age_PPI[:, 0] = np.asarray(list(df['AgeBin']))
    # Age_PPI[:, 1] = np.asarray([df['MedConf_PPI'][i] + df['HighConf_PPI'][i] if df[trait][i] == 1 else 0 for i in df.index])
    Age_PPI[:, 1] = np.asarray([df['HighConf_PPI'][i] if df[trait][i] == 1 else 0 for i in df.index])
    Age_PPI_sorted = Age_PPI[np.argsort(Age_PPI[:, 0])]
    Age_PPI_sorted[np.argwhere(Age_PPI_sorted[:, 1] == 0), 1] = np.nan
    return Age_PPI_sorted


def Allfuncviolins(df):
    Age_BP = np.zeros([len(df), 2])
    Age_BP[:, 0] = np.asarray(list(df['AgeBin']))
    Age_BP[:, 1] = np.asarray([df['MedConf_PPI'][i] + df['HighConf_PPI'][i] for i in df.index])
    # Age_BP[:, 1] = np.asarray([df['HighConf_PPI'][i] for i in df.index])
    Age_BP_sorted = Age_BP[np.argsort(Age_BP[:, 0])]
    Age_BP_sorted[np.argwhere(Age_BP_sorted[:, 1] == 0), 1] = np.nan
    tmp2 = logviolins(Age_BP_sorted, int(max(Age_BP_sorted[:, 0]) + 1))
    return tmp2


def anova_2w(x):
    # get data into single array, pad missing values
    traitCats = ['A', 'M', 'C', 'D', 'I']
    ageCats = ['old', 'young']
    PPIs = []
    Ages = []
    traits = []
    for i in [0, -1]:
        for j in range(len(x[i])):
            numVals = len(x[i][j])
            PPIs = PPIs + x[i][j]
            Ages = Ages + [ageCats[i]]*numVals
            traits = traits + [traitCats[j]]*numVals

    df = pd.DataFrame(columns=['PPIs', 'Ages', 'traits'], data=zip(PPIs, Ages, traits))

    lm = ols('PPIs~C(Ages)+C(traits)', data=df).fit()
    table = sm.stats.anova_lm(lm)
    Age_pval = table['PR(>F)'][['C(Ages)', 'C(traits)']].values

    return Age_pval


def bootstrap(x, steps):
    # find the smallest set of data in x
    smallest = min([len(i) for i in x])
    vals = np.zeros([len(x), steps])
    for i in range(len(x)):
        for j in range(steps):
            vals[i, j] = np.mean(random.sample(x[i], smallest))
    return vals


def func_sep_bins(df, trait):
    tmp = AgeSortedArray(df, trait)
    tmp2 = logviolins(tmp, int(max(tmp[:, 0]) + 1))
    return tmp2


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

# count processes separated by function
Hall = Allfuncviolins(Hdf)
Hmetab = func_sep_bins(Hdf, 'metabolic process')
Hcp = func_sep_bins(Hdf, 'cellular process')
Hdev = func_sep_bins(Hdf, 'developmental process')
Hisp = func_sep_bins(Hdf, 'immune system process')
Hisp = func_sep_bins(Hdf, 'immune system process')
Hdata = list(zip(Hall, Hmetab, Hcp, Hdev, Hisp))
Hpval = anova_2w(Hdata)

Mall = Allfuncviolins(Mdf)
Mmetab = func_sep_bins(Mdf, 'metabolic process')
Mcp = func_sep_bins(Mdf, 'cellular process')
Mdev = func_sep_bins(Mdf, 'developmental process')
Misp = func_sep_bins(Mdf, 'immune system process')
Misp = func_sep_bins(Mdf, 'immune system process')
Mdata = list(zip(Mall, Mmetab, Mcp, Mdev, Misp))
Mpval = anova_2w(Mdata)

Dall = Allfuncviolins(Ddf)
Dmetab = func_sep_bins(Ddf, 'metabolic process')
Dcp = func_sep_bins(Ddf, 'cellular process')
Ddev = func_sep_bins(Ddf, 'developmental process')
Disp = func_sep_bins(Ddf, 'immune system process')
Ddata = list(zip(Dall, Dmetab, Dcp, Ddev))
Dpval = anova_2w(Ddata)

Call = Allfuncviolins(Cdf)
Cmetab = func_sep_bins(Cdf, 'metabolic process')
Ccp = func_sep_bins(Cdf, 'cellular process')
Cdev = func_sep_bins(Cdf, 'developmental process')
Cisp = func_sep_bins(Cdf, 'immune system process')
Cdata = list(zip(Call, Cmetab, Ccp, Cdev))
Cpval = anova_2w(Cdata)

Aall = Allfuncviolins(Adf)
Ametab = func_sep_bins(Adf, 'metabolic process')
Acp = func_sep_bins(Adf, 'cellular process')
Adev = func_sep_bins(Adf, 'developmental process')
Aisp = func_sep_bins(Adf, 'immune system process')
Adata = list(zip(Aall, Ametab, Acp, Adev))
Apval = anova_2w(Adata)

Drall = Allfuncviolins(Drdf)
Drmetab = func_sep_bins(Drdf, 'metabolic process')
Drcp = func_sep_bins(Drdf, 'cellular process')
Drdev = func_sep_bins(Drdf, 'developmental process')
Drisp = func_sep_bins(Drdf, 'immune system process')
Drdata = list(zip(Drall, Drmetab, Drcp, Drdev))
Drpval = anova_2w(Drdata)

pvals = [Hpval, Mpval, Drpval, Dpval, Cpval, Apval]
Data = [Hdata, Mdata, Drdata, Ddata, Cdata, Adata]
colors = ['Purple', 'Blue', 'Orange', 'Red', 'Black']
labels = ['All', 'Met.P', 'Cell.P', 'Dev.P', 'Imm.P']
lineXloc = [4.5, 4.5, 3.5, 3.5, 3.5, 3.5]

pos = []
pos2 = []
for i in range(2):
    pos.append([0, 1, 2, 3, 4])
    pos2.append([5, 6, 7, 8, 9])
for i in range(4):
    pos.append([0, 1, 2, 3])
    pos2.append([4, 5, 6, 7])
ticks = [[2, 7], [2, 7], [1.5, 5.5], [1.5, 5.5], [1.5, 5.5], [1.5, 5.5]]
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(17, 8))
cntr = 0
for ax in axes.flat:
    vp = ax.violinplot(Data[cntr][0], pos[cntr], showmeans=False, showmedians=False, showextrema=False)
    vp2 = ax.violinplot(Data[cntr][-1], pos2[cntr], showmeans=False, showmedians=False, showextrema=False)
    ax.set_xticks(ticks[cntr])
    ax.set_xticklabels(['Oldest', 'Youngest'])
    ax.set_ylim([0, 9])
    ax.plot([lineXloc[cntr], lineXloc[cntr]], [0, 9], color='black')

    cntr2 = 0
    for pc in vp['bodies']:
        pc.set_facecolor(colors[cntr2])
        cntr2 = cntr2 + 1
    cntr2 = 0
    for pc in vp2['bodies']:
        pc.set_facecolor(colors[cntr2])
        cntr2 = cntr2 + 1

    cntr2 = 0
    for i in Data[cntr][0]:
        quartile1, medians, quartile3 = np.percentile(i, [25, 50, 75])
        whiskers = np.array(adjacent_values(i, quartile1, quartile3))
        whiskersMin, whiskersMax = whiskers[0], whiskers[1]

        ax.scatter(pos[cntr][cntr2], medians, marker='_', color='white', s=5, zorder=3)
        ax.vlines(pos[cntr][cntr2], quartile1, quartile3, color='k', linestyle='-', lw=2.5)
        ax.vlines(pos[cntr][cntr2], whiskersMin, whiskersMax, color='k', linestyle='-', lw=.5)
        cntr2 = cntr2 + 1
    cntr2 = 0
    for i in Data[cntr][-1]:
        quartile1, medians, quartile3 = np.percentile(i, [25, 50, 75])
        whiskers = np.array(adjacent_values(i, quartile1, quartile3))
        whiskersMin, whiskersMax = whiskers[0], whiskers[1]

        ax.scatter(pos2[cntr][cntr2], medians, marker='_', color='white', s=5, zorder=3)
        ax.vlines(pos2[cntr][cntr2], quartile1, quartile3, color='k', linestyle='-', lw=2.5)
        ax.vlines(pos2[cntr][cntr2], whiskersMin, whiskersMax, color='k', linestyle='-', lw=.5)
        cntr2 = cntr2 + 1
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
all = mpatches.Patch(color='purple', alpha=.5, label='All')
met = mpatches.Patch(color='blue', alpha=.5, label='Metabolic')
cell = mpatches.Patch(color='orange', alpha=.5, label='Cellular')
Dev = mpatches.Patch(color='red', alpha=.5, label='Developmental')
Imm = mpatches.Patch(color='black', alpha=.5, label='Immune')
fig.legend(handles=[all, met, cell, Dev, Imm], loc='right')
fig.show()
fig.savefig(image_folder+"S2_Age_vs_PPI_boxplot.tiff")


# boot strap means by sub-sampling the values from each functional group by the lowest n group
steps = 50
Hold = bootstrap(Hdata[0], steps)
Hyoung = bootstrap(Hdata[-1], steps)

Mold = bootstrap(Mdata[0], steps)
Myoung = bootstrap(Mdata[-1], steps)

Drold = bootstrap(Drdata[0], steps)
Dryoung = bootstrap(Drdata[-1], steps)

Dold = bootstrap(Ddata[0], steps)
Dyoung = bootstrap(Ddata[-1], steps)

Cold = bootstrap(Cdata[0], steps)
Cyoung = bootstrap(Cdata[-1], steps)

Aold = bootstrap(Adata[0], steps)
Ayoung = bootstrap(Adata[-1], steps)

oldData = [Hold, Mold, Drold, Dold, Cold, Aold]
youngData = [Hyoung, Myoung, Dryoung, Dyoung, Cyoung, Ayoung]
colors = ['Purple', 'Blue', 'Orange', 'Red', 'Black']
labels = ['All', 'Met.P', 'Cell.P', 'Dev.P', 'Imm.P']
pos = []
pos2 = []
for i in range(2):
    pos.append([0, 1, 2, 3, 4])
    pos2.append([5, 6, 7, 8, 9])
for i in range(4):
    pos.append([0, 1, 2, 3])
    pos2.append([4, 5, 6, 7])
ticks = [[2, 7], [2, 7], [1.5, 5.5], [1.5, 5.5], [1.5, 5.5], [1.5, 5.5]]
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(17, 8))
cntr = 0
for ax in axes.flat:
    lw = 2
    cntr2 = 0
    ax.set_xticks(ticks[cntr])
    ax.set_xticklabels(['Oldest', 'Youngest'])
    ax.set_ylim([1, 7])
    ax.plot([lineXloc[cntr], lineXloc[cntr]], [0, 7], color='black')

    for i in oldData[cntr]:
        median = np.mean(i)
        conf = 1.96 * (np.std(i) / np.sqrt(len(i)))
        ax.vlines(pos[cntr][cntr2], median-conf, median+conf, color=colors[cntr2], label=labels[cntr2], linestyle='-', lw=lw)
        ax.hlines(median, pos[cntr][cntr2] - .1, pos[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        ax.hlines(median-conf, pos[cntr][cntr2] - .1, pos[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        ax.hlines(median+conf, pos[cntr][cntr2] - .1, pos[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        cntr2 = cntr2 + 1
    cntr2 = 0
    for i in youngData[cntr]:
        median = np.mean(i)
        conf = 1.96*(np.std(i)/np.sqrt(len(i)))
        ax.vlines(pos2[cntr][cntr2], median-conf, median+conf, color=colors[cntr2], linestyle='-', lw=2)
        ax.hlines(median, pos2[cntr][cntr2] - .1, pos2[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        ax.hlines(median-conf, pos2[cntr][cntr2] - .1, pos2[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        ax.hlines(median+conf, pos2[cntr][cntr2] - .1, pos2[cntr][cntr2] + .1, color=colors[cntr2], linestyle='-', lw=lw)
        cntr2 = cntr2 + 1
    cntr = cntr + 1

axes[0, 0].set_title('H. sapiens')
axes[0, 1].set_title('M. musculus')
axes[0, 2].set_title('D. rerio')
axes[1, 0].set_title('D. melanogaster')
axes[1, 1].set_title('C. elegans')
axes[1, 2].set_title('A. thaliana')
plt.setp(axes[-1, :], xlabel='Age Bin')
plt.setp(axes[0, 0], ylabel='log(PPI)')
plt.setp(axes[1, 0], ylabel='log(PPI)')
all = mpatches.Patch(color='purple', alpha=.5, label='All')
met = mpatches.Patch(color='blue', alpha=.5, label='Metabolic')
cell = mpatches.Patch(color='orange', alpha=.5, label='Cellular')
Dev = mpatches.Patch(color='red', alpha=.5, label='Developmental')
Imm = mpatches.Patch(color='black', alpha=.5, label='Immune')
fig.legend(handles=[all, met, cell, Dev, Imm], loc='right')
fig.show()
fig.savefig(image_folder+"S5_Age_vs_PPI_bootstrap_means.tiff")
