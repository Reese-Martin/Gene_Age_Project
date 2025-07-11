import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
image_folder = 'Images Combined/'
# load Gene age DF
Hdata_folder = 'Data Human/'
Hdf = pd.read_pickle(Hdata_folder+"Gene_Age_DataFrame.pkl")

# count the number of entries there are for each age group
Art_Num_Groups = Hdf['article_count'].to_list()
c = Counter(Art_Num_Groups)

# bar plot of the number of genes there are for each age group
f1 = plt.figure()
plt.bar(c.keys(), c.values())
plt.xlim(0, 1500)
plt.ylim(0, 500)
plt.ylabel('Number of Proteins')
plt.xlabel('Number of Articles')
plt.show()
f1.savefig(image_folder+"Gene_Article_Count.png")
