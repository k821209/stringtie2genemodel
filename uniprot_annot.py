import pandas as pd 
import sys
sys.path.append('/ref/analysis/pipelines')
import kang
# new protein annot with uniprot 
file_bp = sys.argv[1] # 'new_gene.fa.bp.na1.uniprot-all.fasta.out7'
df_bp   = pd.read_csv(file_bp,sep='\t',header=None,comment='#')

file_uniprot = sys.argv[2] #'/ref/analysis/References/uniprot/uniprot-all.fasta'
dic_uniprot  = kang.Fasta2dic_all(file_uniprot)

keys = [x.split()[0] for x in dic_uniprot.keys()]
values = [' '.join(x.split()[1:]) for x in dic_uniprot.keys()]
dic_uniprot_IDmap = dict(zip(keys,values))

df_bp['annotation'] = df_bp[1].apply(lambda x : dic_uniprot_IDmap[x])

df_bp[df_bp[10] < 0.1].to_csv(file_bp+'.annot',sep='\t')

