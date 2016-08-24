from __future__ import print_function
import pandas as pd
import numpy as np
import sys
sys.path.append('/ref/analysis/pipelines/')
import kang
from tqdm import tqdm
import glob

file_stringtie_fa = sys.argv[1] #'/ref/analysis/Cre/braker/braker.try5_mario/guided/transcripts.fasta'
dic_stringtie_fa = kang.Fasta2dic(file_stringtie_fa)

#main_dir      = './' 
file_ag       = sys.argv[2:6] #main_dir+'transcripts.fasta.augustus.ath.complete.gff3.nosharp.genome.v1.gff'


file_td       = sys.argv[6] #main_dir+'selected_mRNA_v4.gff'
ag_predictions = []
for efile_ag in file_ag:
    df_ag     = pd.read_csv(efile_ag,sep='\t')
    df_ag['ID'] = df_ag['Name'].apply(lambda x : x.replace('Name=',''))
    df_ag.set_index('ID',inplace=True)
    ag_predictions.append(df_ag)


df_td     = pd.read_csv(file_td,sep='\t')
#df_td['ID'] = df_td['Name'].apply(lambda x : x.split('|')[0].replace('Name=',''))
df_td['ID'] = df_td['Name'].apply(lambda x : x.split('::')[1])

df_td.set_index('ID',inplace=True)

i = 0 
with open('my_csv.csv', 'w') as f:
    print('Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes','Name',sep='\t',file=f)
    for tid in tqdm(set(dic_stringtie_fa.keys())):
        b = 1
        for edf in ag_predictions:
            try:
                sub_df = edf.loc[tid]
                break
            except KeyError:            
                b = 0 
                pass 


        if b == 1:
            sub_df[['Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes','Name']].to_csv(f,header=False,sep='\t') 
        
        try:
            sub_df = df_td.loc[tid]
            sub_df[['Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes','Name']].to_csv(f,header=False,sep='\t')
        except KeyError:
            continue

