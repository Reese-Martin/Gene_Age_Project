# this file converts the very specific GO terms into broader categories (ie. metabolism)

import pandas as pd
import numpy as np
import networkx
import obonet
from collections import Counter

data_folder = 'Data Dreri/'
common_folder = 'Common Data/'
df = pd.read_pickle(data_folder + 'Gene_Age_DataFrame.pkl')

# load GO terms to a relational graph. The go-basic.obo file was attained from the Gene Ontology website on 11/27/2023
graph = obonet.read_obo(common_folder + 'go-basic.obo')
# generate dictionary of GO ids with name as key
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
Broad_Procs = []
fails = []
for i in df['Assc_BioProcesses']:
    Broad_Proc = []
    for j in i:
        # build paths between the given term and the biological_process
        try:
            paths = networkx.all_simple_paths(graph,
                                              source=name_to_id[j],
                                              target=name_to_id['biological_process']
                                              )
            tmp = []
            for path in paths:
                tmp.append(''.join(id_to_name[node] for node in path[-2:-1]))
            # return the most prominent process from those collected
            counts = Counter(tmp)
            Broad_Proc.append(max(counts, key=counts.get))
        except:
            fails.append(j)
    Broad_Procs.append(Broad_Proc)
df['Broad_BioProcesses'] = Broad_Procs
df.to_pickle(data_folder+'Gene_Age_DataFrame.pkl')
