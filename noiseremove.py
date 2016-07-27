import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import sys
matplotlib.style.use('ggplot') 

file_in = sys.argv[1] # 'transcripts.fasta.transdecoder.gff3.genome.gff'
df      = pd.read_csv(file_in,sep='\t')
df['transcript_name'] = df['GFF3: grouping attributes'].apply(lambda x : x.split(';')[-1].replace('Name=',''))

# ratio of UTR length to total length of ORF 
mask   = (df['feature'] == 'five_prime_UTR' ) | (df['feature'] == 'three_prime_UTR' ) | (df['feature'] == 'CDS' )
df_sub = df[mask]
df_sub['transcript_name'] = df_sub['GFF3: grouping attributes'].apply(lambda x : x.split(';')[-1].replace('Name=',''))
df_sub_ix                 = df_sub.set_index('transcript_name')

def covering(mat):
    cov = 0
    for a,b in np.array(mat):
        cov += b-a
    return cov
def covering_g(mat):
    if len(mat) == 0 :
        return 0 
    else:
        return np.max(mat) - np.min(mat)


dic = {}
for i in list(set(df_sub_ix.index)):
    edf = df_sub_ix.loc[i]
    if isinstance(edf,pd.Series):
        if edf['feature'] != 'CDS':
            print edf
            break
        left,right = edf[['start','end']].values
        dic[i] = 100,100
    else:
        cds_locs = edf[edf['feature'] == 'CDS'][['start','end']].values
        utr3_locs = edf[(edf['feature'] == 'three_prime_UTR')][['start','end']].values
        utr5_locs = edf[(edf['feature'] == 'five_prime_UTR')][['start','end']].values
        
        #cds_cov, utr_cov = covering(cds_locs),covering(utr_locs)
        cds_cov, utr_cov = covering_g(cds_locs),covering_g(utr3_locs)+covering_g(utr5_locs)
        dic[i] = float(utr_cov)/(cds_cov+utr_cov),cds_cov,


o = np.array([x[0] for x in dic.values()])
transcript_list = np.array(dic.keys())[o<0.4]
def is_member(a):
    return a in transcript_list
is_member_vec = np.vectorize(is_member)
mask = is_member_vec(df['transcript_name'])

selected_df = df[mask]
selected_df.to_csv('selected_mRNA_v4.gff',sep='\t')

