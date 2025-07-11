
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def average_count(x):
    tmp = list(set(x[:, 0]))
    tmp = [i for i in tmp if ~np.isnan(i)]
    tmp.sort()
    array = np.zeros([len(tmp), 2])
    for i in range(len(tmp)):
        array[i, 0] = tmp[i]
        array[i, 1] = np.mean(x[np.argwhere(x[:, 0] == tmp[i]), 1])
    return array


def windowed_average(time_series, window_size):
    result = []
    n = len(time_series)

    for i in range(n):
        # Determine the window range
        start = max(0, i - window_size + 1)  # Ensure we don't go negative
        end = i + 1

        # Calculate the windowed average
        window = time_series[start:end]
        if len(window) < window_size:
            # If the window is smaller than expected, return the raw value
            result.append(time_series[i])
        else:
            # Calculate the average of the current window
            average = sum(window) / len(window)
            result.append(average)

    return result


image_folder = 'Images Combined/'
# load Gene age DF
Hdata_folder = 'Data Human/'
Hdf = pd.read_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")

# article count vs BioProcesses
AC_BP = np.zeros([len(Hdf), 2])
AC_BP[:, 0] = np.asarray(list(Hdf['article_count']))
AC_BP[:, 1] = np.asarray([len(Hdf.loc[i,'Assc_BioProcesses']) for i in Hdf.index])
data = average_count(AC_BP)
data[:, 1] = windowed_average(data[:, 1], 40)

f1 = plt.figure()
plt.scatter(AC_BP[:, 0], AC_BP[:, 1])
plt.plot(data[:, 0], data[:, 1], c='orange')
plt.ylabel('Number of Biological Processes')
plt.xlabel('number of Articles')
plt.show()
f1.savefig(image_folder+"Article_Count_vs_BP_Count.png")

AC_PPI = np.zeros([len(Hdf), 2])
AC_PPI[:, 0] = np.asarray(list(Hdf['article_count']))
# enter low confidence counts as the y-value
AC_PPI[:, 1] = Hdf['UHC_PPI'].values
data = average_count(AC_PPI)
data[:, 1] = windowed_average(data[:, 1], 40)

f1 = plt.figure()
plt.scatter(AC_PPI[:, 0], AC_PPI[:, 1])
plt.plot(data[:, 0], data[:, 1], c="orange")
plt.ylabel('Number of Ultra High Confidence PPI')
plt.xlabel('Number of articles')
plt.show()
f1.savefig(image_folder+"Article_Count_vs_UHC_PPI_Count.png")
