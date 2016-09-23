
from __future__ import print_function
import pandas as pd
import numpy as np
import sys
sys.path.append('/ref/analysis/pipelines/')
import kang


file_tmap     = sys.argv[1] 
file_gff3     = sys.argv[2] 
file_gff3_old = sys.argv[3] 
file_pep_new  = sys.argv[4] 

df_tmap   = pd.read_csv(file_tmap,sep='\t')
df_gff = pd.read_csv(file_gff3,sep='\t',header=None,comment='#')
df_gff_old = pd.read_csv(file_gff3_old,sep='\t',header=None,comment='#')

def get_name(feature,x):
    key   = [i.split('=')[0] for i in x.split(';')]
    value = [i.split('=')[1] for i in x.split(';')]
    dic   = dict(zip(key,value))
    return dic[feature]


mask       = (df_gff[2] == 'mRNA')
df_gff_mRNA = df_gff[mask]
df_gff_mRNA['Name'] = df_gff_mRNA[8].apply(lambda x : get_name('Name',x))

mask       = (df_gff[2] != 'mRNA') & (df_gff[2] != 'gene')
df_gff_CDS = df_gff[mask]
df_gff_CDS['Name'] = df_gff_CDS[8].apply(lambda x : get_name('Name',x))

df_gff_transcript = df_gff_mRNA.append(df_gff_CDS)
df_gff_transcript_ix = df_gff_transcript.set_index('Name')

mask        = (df_gff[2] == 'gene') 
df_gff_gene = df_gff[mask]
df_gff_gene['Name'] = df_gff_gene[8].apply(lambda x : get_name('ID',x))

mask        = (df_gff[2] == 'mRNA') 
df_gff_mRNA = df_gff[mask]
df_gff_mRNA['Name'] = df_gff_mRNA[8].apply(lambda x : get_name('Parent',x))



mask        = (df_gff_old[2] == 'mRNA')
df_gff_old_mRNA = df_gff_old[mask]
df_gff_old_mRNA['Name'] = df_gff_old_mRNA[8].apply(lambda x : get_name('Parent',x))
df_gff_old_mRNA['tName'] = df_gff_old_mRNA[8].apply(lambda x : get_name('ID',x))
dicT2G      = dict(zip(df_gff_old_mRNA['tName'], df_gff_old_mRNA['Name']))






mask       = (df_gff[2] == 'CDS')
df_gff_CDS = df_gff[mask]
df_gff_CDS['Name'] = df_gff_CDS[8].apply(lambda x : '.'.join(get_name('Name',x).split('.')[0:2]))

df_gff_genename = df_gff_gene.append(df_gff_mRNA.append(df_gff_CDS))
df_gff_genename_ix = df_gff_genename.set_index('Name')

# adding isoforms 
mask           = (df_tmap['class_code'] == 'j')
df_tmap_sub    = df_tmap[mask]
df_tmap_sub_ix = df_tmap_sub.set_index('cuff_id')

with open(file_gff3+'.merge.gff3','w') as f:
    for ix in df_tmap_sub_ix.index:
        preparent = '.'.join(ix.split('.')[0:2])
        info      = 'ID=%s;Parent=%s;Name=%s'%(ix,dicT2G[df_tmap_sub_ix.loc[ix]['ref_id']],ix)
        edf       =  df_gff_transcript_ix.loc[ix]
        mRNA_df   = edf[edf[2] == 'mRNA']
        Other_df  = edf[edf[2] != 'mRNA']
        print ('\t'.join(map(str,mRNA_df[[0,1,2,3,4,5,6,7]].values[0])),info,sep='\t',file=f)
        Other_df.to_csv(f,header=None,index=None,sep='\t')  

# adding new genes 
mask           = (df_tmap['class_code'] == 'u') | (df_tmap['class_code'] == 'o')
df_tmap_sub    = df_tmap[mask]
df_tmap_sub_ix = df_tmap_sub.set_index('cuff_id')

genelist = []

with open(file_gff3+'.merge.gff3','a') as f:
    for ix in set(df_tmap_sub_ix.index):
        genename = '.'.join(ix.split('.')[0:2])
        if genename not in set(genelist):
            edf       = df_gff_genename_ix.loc[genename]
            edf       = edf[edf[2] =='gene']
            edf.to_csv(f,header=None,index=None,sep='\t')    
            genelist.append(genename)
        else:
            pass
        edf       = df_gff_transcript_ix.loc[ix]
        edf.to_csv(f,header=None,index=None,sep='\t')    

# new gene protein seq ret 
file_fa       = file_pep_new
dicfa         = kang.Fasta2dic(file_fa)

df_new        = df_tmap_sub
df_new['seq'] = df_new['cuff_id'].apply(lambda x : dicfa[x])

with open(file_pep_new+'.new_gene.fa','w') as f:
    for ix in set(df_new.index):
        hd  = df_new.loc[ix][11]
        seq = df_new.loc[ix]['seq']
        print('>'+hd,file=f)
        print(seq,file=f)
# new gene protein seq ret end 
