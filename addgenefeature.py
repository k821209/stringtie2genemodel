#!/usr/bin/python
'''
        Chromosome      source  feature start   end     score   strand  frame   GFF3: grouping attributes       Name
288813  chromosome_1    ST/TD/KYJ       mRNA    245102  248187  .       -       .       ID=STRG.10.1|m.23339;Parent=STRG.10.1;Name=STRG.10.1|m.23339    Name=STRG.10.1|m.23339
288816  chromosome_1    ST/TD/KYJ       five_prime_UTR  245102  245891  .       -       .       ID=STRG.10.1|m.23339.utr5p1.0;Parent=STRG.10.1|m.23339;Name=STRG.10.1|m.23339   Name=STRG.10.1|m.23339
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import sys
file_in = sys.argv[1]#'test.augustus.gff3.nosharp.gff3.genome.gff' 
df      = pd.read_csv(file_in,sep='\t')
mask    = (df['feature'] == 'mRNA')
sub_df  = df[mask]

sub_df['geneid'] = sub_df['Name'].apply(lambda x : '.'.join(x.replace('Name=','').split('.')[0:2]))
#sub_df['geneid'] = sub_df['transcript_name'].apply(lambda x : '.'.join(x.replace('Name=','').split('.')[0:2]))
#sub_df['geneid'] = sub_df['GFF3: grouping attributes'].apply(lambda x : '.'.join(x.split(';')[-1].replace('Name=','').split('.')[0:2]))
sub_df.set_index('geneid',inplace=True)

adddic = {
         'Chromosome' : [] ,
         'source'     : [] , 
         'feature'    : [] , 
         'start'      : [] , 
         'end'        : [] ,
         'score'      : [] ,
         'strand'     : [] ,
         'frame'      : [] , 
         'GFF3: grouping attributes'  : [] , 
         'Name' : []     
         }
for ix in set(sub_df.index):
    mat        = sub_df.loc[ix][['start','end']].values
    try:
        chromosome = sub_df.loc[ix]['Chromosome'].values[0]
        strand     = sub_df.loc[ix]['strand'].values[0]
    except AttributeError:
        chromosome = sub_df.loc[ix]['Chromosome']
        strand     = sub_df.loc[ix]['strand']
    left       = np.min(mat)
    right      = np.max(mat)
    info       = 'ID='+ix
    name       = 'Name='+ix
    adddic['Chromosome'].append(chromosome)
    adddic['source'].append('ST/ALL/KYJ')
    adddic['feature'].append('gene')
    adddic['start'].append(left)
    adddic['end'].append(right)
    adddic['score'].append('.')
    adddic['strand'].append(strand)
    adddic['frame'].append('.')
    adddic['GFF3: grouping attributes'].append(info)
    adddic['Name'].append(name) 
add_df = pd.DataFrame(adddic,columns=['Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes','Name'])
df_new = df.append(add_df,ignore_index=True)

Outfile = open(file_in+'.addgene.gff3','w')
print('##gff-version 3',file=Outfile)
df_new[['Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes']].sort_values(by=['Chromosome','start']).to_csv(Outfile,sep='\t',header=None,index=False)
#add_df.to_csv(file_in+'.gene.gff',sep='\t')

    
    
    
    



