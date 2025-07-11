import pandas as pd
import os
os.environ['NCBI_API_KEY'] =  # removed for publication as the NCBI API key is private
from metapub import PubMedFetcher
fetch = PubMedFetcher()

df = pd.read_pickle('Gene_Age_DataFrame.pkl')
Prot_Names = df.Full_Name.to_list()

Num_Articles = []
results_list = pd.DataFrame(columns=['Results_count'], index=df.index)

for i in df.index[19918:]:
    if i % 100 == 0:
        print(i/len(Prot_Names))
    if type(Prot_Names[i]) == str:
        results_list['Results_count'][i] = len(fetch.pmids_for_query(Prot_Names[i], retmax=1500))

df['article_count'] = results_list['Results_count']
df.to_pickle('Gene_Age_DataFrame.pkl')

