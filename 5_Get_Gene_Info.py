
import pandas as pd
from omadb import Client
import stringdb
import numpy as np

data_folder = 'Data Mouse/'

# load Gene age DF
df = pd.read_pickle(data_folder+"Gene_Age_DataFrame.pkl")

# set up OMAdb client
c = Client()

# lists that store the number of biological processes each gene is associated with
# there are a ton of evidence classes, so for now we will just use all evidence
# this can take a long time (> 2 hours for ~20k genes)
All_evidence = []  # All processes
Exp_evidence = []
Missed = []
timeouts = []  # record when the API timed out in case we need to re-trace steps
cntr = 0
evidences = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']

for i in df['Gene'][15335:]:
    tmp = c.entries.gene_ontology(i, aspect="BP")
    if tmp:  # checks to make sure there was a biological process associated with the given protein
        tmp2 = []
        for j in tmp[i[0]]:
            if j['evidence'] in evidences:
                tmp2.append(j['name'])
        Exp_evidence.append(list(set(tmp2)))
        All_evidence.append(list(set([k['name'] for k in tmp[i[0]]])))
    else:
        Exp_evidence.append(tmp)
        All_evidence.append(tmp)
        Missed.append(cntr)
    cntr = cntr+1
    if cntr % 100 == 0:
        print(cntr/len(df))

### testing
test = df['Assc_BioProcesses'].tolist()
val = 19664
print(f"original {test[val]} and new {All_evidence[val]}")

df['Assc_BioProcesses'] = All_evidence
df['Assc_BioProcesses_Exp'] = Exp_evidence

# generate string PPI for each gene
stringDB_names = []
timeouts = []
cntr = 0
# due to API timeouts on the part of OMADB, you may need to run this multiple times
# updating the start point to be the last entry from the timeout
for i in df['Gene']:
    try:
        stringDB_names.append(c.entries[i][0]['canonicalid'])
    except:
        print(' API timeout')
        timeouts.append(i)
        stringDB_names.append(['TO'])
    cntr = cntr+1
    if cntr % 100 == 0:
        print(cntr/len(df))
df['String_Name'] = stringDB_names

indices = [i for i in df.index.to_list() if df['String_Name'][i] == ['TO']]  # captures time out indices
for i in range(0, len(indices)):
    tmp = c.entries[timeouts[i]]
    if tmp:
        df.loc[indices[i], 'String_Name'] = tmp[0]['canonicalid']
    else:
        df.at[indices[i], 'String_Name'] = []

# QC to remove the genes that have no string name
indices = [i for i in df.index.to_list() if df['String_Name'][i] == []]  # captures time out indices
df.drop(indices, inplace=True)
df = df.reset_index(drop=True)
# retrieve string IDs to be used for searches in the interaction dataframes
# API timeout if you request too many more genes than 10k so break queries into chunks
# species codes:
# Human: 9606
# Celeg: 6239
# Mouse: 10090
# Dmel: 7227
# Athal: 3702
# Dreri: 7955
species = 3702
String_IDs_1 = stringdb.get_string_ids(df.String_Name[:10000], species=species)
String_IDs_2 = stringdb.get_string_ids(df.String_Name[10000:20000], species=species)
# this last line is only necessary for the mouse dataset, Athaliana, and d rerio
String_IDs_3 = stringdb.get_string_ids(df.String_Name[20000:], species=species)
String_IDs = pd.concat([String_IDs_1,  String_IDs_2, String_IDs_3], ignore_index=True)
# String_IDs = pd.concat([String_IDs_1,  String_IDs_2], ignore_index=True)
df['String_IDs'] = ''
# assign collected string IDs to the associated protein, cannot just set the arrays equal, because some
# entries do not have a String ID
for i in range(0, len(String_IDs)):
    df.loc[df['String_Name'] == String_IDs.queryItem[i], 'String_IDs'] = String_IDs.loc[i, 'stringId']

# using string IDs, count the number of PPI from the low, medium, and high PPI dataframes
#  this takes a while,

data_folder = 'Data Mouse/'

# load Gene age DF
df = pd.read_pickle(data_folder+"Gene_Age_DataFrame.pkl")
df['UHC_PPI'] = np.zeros(len(df))
# df["MedConf_PPI"] = ''
# df["HighConf_PPI"] = ''
# MedConf_df = pd.read_pickle(data_folder+"Med_Confidence_string_Interactions.pkl")
# HighConf_df = pd.read_pickle(data_folder+"High_Confidence_string_Interactions.pkl")
UHC_df = pd.read_pickle(data_folder+"Ultra_High_Confidence_string_Interactions.pkl")
# this can be slow, but you should only need to run it once. Just be sure to save the data

for i in range(0, len(df)):
    if i % 100 == 0:
        print(i/len(df))
    tmp = df['String_IDs'][i]
    # df.loc[i, 'MedConf_PPI'] = len(MedConf_df[MedConf_df['protein1'] == tmp])
    # df.loc[i, 'HighConf_PPI'] = len(HighConf_df[HighConf_df['protein1'] == tmp])
    df.loc[i, 'UHC_PPI'] = len(UHC_df[UHC_df['protein1'] == tmp])

# This region assigns human readable gene symbols to each gene based on the string ID of a protein
# Human: 9606
# Mouse: 10090
# Dreri: 7955
# Dmel: 7227
# Celeg: 6239
# Athal: 3702
Spec = '7955'
Full_Name_map = pd.read_csv(data_folder+Spec+".protein.info.v11.5.txt", sep='\t')
Full_GeneNames = []
for i in range(0, len(df)):
    tmp = Full_Name_map[Full_Name_map['#string_protein_id'] == df['String_IDs'][i]]
    if tmp['preferred_name'].to_list():
        Full_GeneNames.append(tmp['preferred_name'].to_list()[0])
    else:
        Full_GeneNames.append('')
df['Gene_Symbol'] = Full_GeneNames
# first get the gene names from uniprot IDs so that we can appropriately access the duplicated genes
# we do this by saving the uniprot id's to a txt file, then uploading this file to the https://www.uniprot.org/id-mapping
# then from ID mapping file we get the gene name
UniprotIDs = df['String_Name'].to_list()
# open file in write mode
with open(data_folder+'UniProtIDs.txt', 'w') as fp:
    for item in UniprotIDs:
        # write each item on a new line
        fp.write("%s\n" % item)
    print('Done')

# using data pulled from ensemble biomart to assign paralogy
# load Gene age DF
# df = pd.read_pickle(data_folder+"Gene_Age_DataFrame.pkl")
df.rename(columns={'Duplication': 'Paralogs'}, inplace=True)
df['Paralogs'] = 0
paralogs = pd.read_csv(data_folder+'mart_export.txt', index_col=0)
paralogs.dropna(inplace=True)

dfgenes = df['Gene_Symbol'].tolist()
counts = np.unique(paralogs['Gene name'], return_counts=True)

for i, j in zip(counts[0], counts[1]):
    if i in dfgenes:
        df.loc[df['Gene_Symbol'] == i, ['Paralogs']] = j

test = np.array(df['Paralogs'].tolist())
np.count_nonzero(test)
df.to_pickle(data_folder+"Gene_Age_DataFrame.pkl")
