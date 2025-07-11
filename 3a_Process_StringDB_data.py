# this script takes the downloaded string data and breaks it apart into the low, medium, and high confidence
# string scores so that string interactions to gene age comparisons can be made independently
import pandas as pd
common_folder = 'Common Data/'
data_folder = 'Data Human/'

# load string txt to dataframe
df = pd.read_csv(data_folder+'4932.protein.links.v12.0.txt', sep=" ")
df['combined_score'] = df['combined_score'].values/1000

med_conf_Dicts = []
high_conf_Dicts = []
for i in range(0, len(df.combined_score.values)):
    if (df.iloc[i, 2] < .7) & (df.iloc[i, 2] >= .4):
        med_conf_Dicts.append({'protein1': df.iloc[i, 0], 'protein2': df.iloc[i, 1], 'combined_score': df.iloc[i, 2]})
    elif df.iloc[i, 2] >= .7:
        high_conf_Dicts.append({'protein1': df.iloc[i, 0], 'protein2': df.iloc[i, 1], 'combined_score': df.iloc[i, 2]})

med_conf = pd.DataFrame.from_dict(med_conf_Dicts)
high_conf = pd.DataFrame.from_dict(high_conf_Dicts)

med_conf.to_pickle(data_folder+"Med_Confidence_string_Interactions.pkl")
high_conf.to_pickle(data_folder+"High_Confidence_string_Interactions.pkl")

# using the high confidence data, extract only the interactions with .95 or greater scores
data_folder = 'Data Celeg/'
df = pd.read_pickle(data_folder+"High_Confidence_string_Interactions.pkl")
print(len(df),'\n')
UHC_df = df[df['combined_score'] >= .95]
print(len(UHC_df))
UHC_df.to_pickle(data_folder+"Ultra_High_Confidence_string_Interactions.pkl")