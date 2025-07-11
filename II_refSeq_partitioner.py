# this script takes in the OMAdb OMA ID to refseq document, and then breaks it up into the individual
# organisms of interest. Saving time when searching refSeq values

import pandas as pd

file = pd.read_csv('Common Data/oma-refseq.txt', sep='\t')

mouse = file[[True if 'MOUSE' in i else False for i in file['OMA ID']]]
mouse.to_csv('Data Mouse/oma-refseq.csv')

human = file[[True if 'HUMAN' in i else False for i in file['OMA ID']]]
human.to_csv('Data Human/oma-refseq.csv')

dmel = file[[True if 'DROME' in i else False for i in file['OMA ID']]]
dmel.to_csv('Data Dmel/oma-refseq.csv')

danio = file[[True if 'DANRE' in i else False for i in file['OMA ID']]]
danio.to_csv('Data Dreri/oma-refseq.csv')

athal = file[[True if 'ARATH' in i else False for i in file['OMA ID']]]
athal.to_csv('Data Athal/oma-refseq.csv')

celeg = file[[True if 'CAEEL' in i else False for i in file['OMA ID']]]
celeg.to_csv('Data Celeg/oma-refseq.csv')
