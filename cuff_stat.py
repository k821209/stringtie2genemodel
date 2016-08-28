#!/usr/bin/python2.7

from __future__ import print_function
import pandas as pd
import sys
sys.path.append('/ref/analysis/pipelines/')
import kang
from collections import Counter 
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

try:
    file_tmap        = sys.argv[1] #'./predicted/cuffcmp.my_csv.csv.addgene.gff3.sort.gff3.tmap'
except IndexError:
    print(''' args : tmap, refcds, pep, cds ''')
    exit()
file_ref_cds     = sys.argv[2] #'/ref/analysis/References/Creinhardtii/annotation/Creinhardtii_281_v5.5.cds.fa'
file_pep         = sys.argv[3] #'../gff2cds/pep.fa'
file_cds         = sys.argv[4] #'../gff2cds/cds.fa'

dicrefcds   = kang.Fasta2dic(file_ref_cds)
dicpep      = kang.Fasta2dic(file_pep)
diccds      = kang.Fasta2dic(file_cds)

df_tmap     = pd.read_csv(file_tmap,sep='\t',comment='#')

mask        = (df_tmap['class_code'] == '=')
TID_ws      = set(df_tmap['ref_id'][mask])  # ws : well supported 
GID_ws      = set(df_tmap['ref_gene_id'][mask])

mask        = (df_tmap['class_code'] == 'c')
TID_cs      = set(df_tmap['ref_id'][mask])  # cs : partial supported
GID_cs      = set(df_tmap['ref_gene_id'][mask])

mask         = (df_tmap['class_code'] == 'j')
TID_js       = set(df_tmap['ref_id'][mask])  # js : isoform supported
GID_js       = set(df_tmap['ref_gene_id'][mask])
TID_js_count = Counter(df_tmap['ref_id'][mask])

mask               = (df_tmap['class_code'] == 'u')
cuff_TID_new       = set(df_tmap['cuff_id'][mask])  # new genes 


with open('cuff_stat.txt','w') as f:
    print('## only perpect support (=), and isoform suport (j) are considered for this summary',file=f)
    print('a. number of genes perpectly  supported            : %d'%len(GID_ws),file=f)
    print('b. number of transcripts perpectly supported       : %d'%len(TID_ws),file=f)
    print('c. number of transcripts only supported as isoform : %d'%len(TID_js-TID_ws-TID_cs),file=f)
    print('d. number of new transcripts                       : %d'%len(cuff_TID_new),file=f)
# generate histogram of number of isoforms on each ref_tid

fig, ax = plt.subplots(1)
ax.hist(TID_js_count.values(),bins=50)
ax.set_xlabel('number of isoforms')
ax.set_ylabel('number of reference transcripts')
#plt.ylim(0,3000)
plt.savefig('hist.isoform.svg')


# generate pie chart

fig, ax = plt.subplots(1,figsize=(10,10))
cuff_class = {'=' : 'Complete match', 'c' : 'Contained (partial match)', 'j' : 'Potentially novel isoform', 'e' : 'possible pre-mRNA fragment', 'i' :  'within a reference intron',
              'o' : 'Generic exonic overlap', 'p' : 'Possible polymerase run-on fragment', 'r' : 'Repeat', 'u': 'Unknown', 'x' : 'opposite strand', 's' : 'opposite strand intron'}

labels = Counter(df_tmap['class_code']).keys()
labels.sort(key=lambda x:Counter(df_tmap['class_code'])[x],reverse=True)
sizes  = Counter(df_tmap['class_code']).values()
sizes.sort(reverse=True)
labels = [cuff_class[x] for x in labels]
n=len(labels) # number of colors
colors=cm.rainbow(np.linspace(0,1,n))
explode = np.zeros(len(labels))
explode[-1] = 0.5
explode[-2] = 0.4
explode[-3] = 0.3
explode[-4] = 0.2
explode[-5] = 0.1
ax.pie(sizes,  
        explode=explode,
        labels=labels,colors=colors,
        autopct='%1.1f%%', shadow=True, startangle=90)
# Set aspect ratio to be equal so that pie is drawn as a circle.
plt.axis('equal')
textstr = '''* Total ref transcripts in gene model : %d, 
  Total ref transcripts analyzed : %d, 
  Total identified transcripts : %d 
  The percentages based on total identified transcripts'''%(len(dicrefcds.keys()),len(set(df_tmap['ref_id'])),len(set(df_tmap['cuff_id'])))
ax.text(-1.6, 1.2, textstr,  fontsize=14,
        verticalalignment='top')
#plt.gcf().subplots_adjust(right=0.05)
plt.gcf().subplots_adjust(left=0.1)
plt.tight_layout()
#plt.gcf().tight_layout()
plt.savefig('pie.chart.tmap.svg',bbox_inches='tight' )







