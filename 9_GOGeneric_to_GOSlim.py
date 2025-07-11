# convert go terms to go slim terms using the obos for both ontologies

import pandas as pd
import obonet
from goatools import obo_parser
from goatools.mapslim import mapslim

# Example Usage
common_folder = 'Common Data/'
data_folder = 'Data Human/'
df = pd.read_pickle(data_folder + 'Gene_Age_DataFrame.pkl')

# load GO terms to a relational graph. The go-basic.obo file was attained from the Gene Ontology website on 11/27/2023
go_dag = obo_parser.GODag(common_folder + 'go-basic.obo')
goslim_dag = obo_parser.GODag(common_folder + 'goslim_generic.obo')

Basic_graph = obonet.read_obo(common_folder + 'go-basic.obo')
Slim_graph = obonet.read_obo(common_folder + 'goslim_generic.obo')
# generate dictionary of GO ids with name as key
Basic_name_to_id = {data['name']: id_ for id_, data in Basic_graph.nodes(data=True) if 'name' in data}

Slim_id_to_name = {id_: data.get('name') for id_, data in Slim_graph.nodes(data=True)}

procs = df['Assc_BioProcesses'].tolist()
slim_procs = []
BasicKeys = Basic_name_to_id.keys()
for proc in procs:
    tmp = []
    for i in proc:
        if i in BasicKeys:
            id_set = mapslim(Basic_name_to_id[i],  go_dag, goslim_dag)[1]
            if id_set:
                tmp.append(Slim_id_to_name[id_set.pop()])
    slim_procs.append(list(set(tmp)))

df['Slim_BioProcesses'] = slim_procs
df.to_pickle(data_folder + 'Gene_Age_DataFrame.pkl')