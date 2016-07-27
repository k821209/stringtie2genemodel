SRRIDs = ["SRR1617532","SRR1617533","SRR1617534","SRR1617535"]
INDEX  = "ref/Esalsugineum_173_v1"
FA     = "ref/Esalsugineum_173_v1.fa"
PFAM   = "/ref/analysis/References/Pfam-A.hmm"
GFF    = "ref/Esalsugineum_173_v1.0.gene.gff3"
UNIPROT = "/ref/analysis/References/uniprot/uniprot-all.fasta"

rule end:
     input : "finalout/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3.sort.gff3","finalout/newgene.annot"
           
# Single end 
rule Hisat2:
     input  : 
             fwd="data/{SRRID}.fastq.gz"
     params : ix=INDEX
     output : 
             "mapped/{SRRID}.bam"
     threads : 10
     shell  : 
             "/program/hisat2-2.0.3-beta/hisat2 --max-intronlen 30000 -p {threads} -x {params.ix} -U {input.fwd}  | sambamba view -f bam -o {output} -S /dev/stdin"

# add "--rna-strandness FR #RF for dUTP protocol" if your library is constructed based on strandness 

rule sambamba_sort : 
     input  : 
            "mapped/{SRRID}.bam"
     output :
            "mapped/{SRRID}.sorted.bam",
            "mapped/{SRRID}.sorted.bam.bai"
     threads : 10
     shell  : 
            "sambamba sort -t {threads} {input}"

#Usage: sambamba-merge [options] <output.bam> <input1.bam> <input2.bam> [...]
rule sambamba_merge : 
     input  : 
            bam=expand("mapped/{SRRID}.sorted.bam",SRRID=SRRIDs),
            bai=expand("mapped/{SRRID}.sorted.bam.bai",SRRID=SRRIDs)
     output : 
            "merged/all.merged.bam"
     threads : 10
     shell  : "sambamba merge -t {threads} {output} {input.bam}"
            
rule stringtie:
     input  :
            "merged/all.merged.bam"
     output : 
            "st_gff_out/all.merged.bam.stringtie.gff"
     threads : 10 
     shell  : 
            "/program/stringtie-1.2.2.Linux_x86_64/stringtie -p {threads} -o {output} {input}"


rule stringtie2gff3:
     input  :
            "st_gff_out/all.merged.bam.stringtie.gff"
     output : 
            "st_gff_out/all.merged.bam.stringtie.gff3"
     shell  : 
            "/program/TransDecoder-2.1.0/util/cufflinks_gtf_to_alignment_gff3.pl {input} > {output}"

rule stringtie2cds:
     input  : "st_gff_out/all.merged.bam.stringtie.gff"
     output : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     params : fa=FA
     shell  : 
            "/program/TransDecoder-2.1.0/util/cufflinks_gtf_genome_to_cdna_fasta.pl {input}  {params.fa} > {output}"


### transdecoder 

rule transdecoder_longOrfs:
     input  : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     output : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/longest_orfs.pep" 
     shell  : 
            "/program/TransDecoder-2.1.0/TransDecoder.LongOrfs -t {input}"

rule hmm:
     input  : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/longest_orfs.pep" 
     output : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/pfam.domtblout"
     threads: 8
     params : pfam=PFAM
     shell  : 
            "hmmscan --cpu {threads} --domtblout {output} {params.pfam} {input}"

rule transdecoder_predict:
     input  : fa="st_gff_out/all.merged.bam.stringtie.gff.fa", pfam="all.merged.bam.stringtie.gff.fa.transdecoder_dir/pfam.domtblout"
     output : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3" 
     threads: 10
     shell  : 
            '''/program/TransDecoder-2.1.0/TransDecoder.Predict --cpu {threads} -t {input.fa} --retain_pfam_hits {input.pfam}        
               mv all.merged.bam.stringtie.gff.fa.transdecoder.* predicted/'''
        
rule transdecoder_togenome:
     input  : sgff="st_gff_out/all.merged.bam.stringtie.gff", tgff="predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3"
     output : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3.genome.gff"
     shell  : "python2.7 stringtie.addcds.py {input.sgff} {input.tgff}"

rule transdecoder_noiseremove:
     input  : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3.genome.gff"
     output : "predicted/selected_mRNA_v4.gff"
     shell  : '''python2.7 noiseremove.py {input} 
                 mv selected_mRNA_v4.gff predicted/'''
###

### Augustus 

rule augustus:
     input  : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     output : "predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3"
     shell  : "/program/augustus-3.2.1/bin/augustus --species=arabidopsis --genemodel=complete --gff3=on --strand=forward {input} > {output}"

rule augustus_togenome:
     input  : sgff="st_gff_out/all.merged.bam.stringtie.gff",tgff="predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3"
     output : "predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff"
     shell  :  '''python2.7 stringtie.augustus.addcds.py {input.sgff} {input.tgff}
                  mv all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff '''

###

rule integrate:
     input  : fa="st_gff_out/all.merged.bam.stringtie.gff.fa", ag="predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff", td="predicted/selected_mRNA_v4.gff"
     output : "predicted/my_csv.csv"
     shell  : '''python2.7 integrate.py {input.fa} {input.ag} {input.td}
                 mv my_csv.csv predicted/'''


rule addgenes:
     input  : "predicted/my_csv.csv"
     output : "predicted/my_csv.csv.addgene.gff3"
     shell  : "python2.7 addgenefeature.py {input}" 
                
rule gt_sort:
     input  : "predicted/my_csv.csv.addgene.gff3"
     output : "predicted/my_csv.csv.addgene.gff3.sort.gff3"
     shell  : "gt gff3 -sort -retainids {input} > {output}" 


rule cuffcompare:
     input  : "predicted/my_csv.csv.addgene.gff3.sort.gff3"
     output : "predicted/cuffcmp.my_csv.csv.addgene.gff3.sort.gff3.tmap"
     params : rgff=GFF
     shell  : '''/program/cufflinks-2.2.1.Linux_x86_64/cuffcompare -r {params.rgff} {input}
                 mv cuffcmp.* predicted/'''

rule merge:
     input  : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3", tmap="predicted/cuffcmp.my_csv.csv.addgene.gff3.sort.gff3.tmap", gff_old=GFF, pepfa_new="finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa"
     output : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3", fa="finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa.new_gene.fa"           
     shell  : '''python2.7 merge.py {input.tmap} {input.gff} {input.gff_old} {input.pepfa_new}
               cat ref/Esalsugineum_173_v1.0.gene.gff3 predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.gff3 > {output.gff}'''

rule gt_sort2:
     input  : "predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3"
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3.sort.gff3"
     shell  : "gt gff3 -sort -retainids {input} > {output}"


rule gff2cds:
     input  : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3"
     params : fa=FA
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.cds.fa"
     shell  : "python2.7 gff2cds.py {params.fa} {input.gff} {output}"

rule cds2pep:
     input  : "finalout/my_csv.csv.addgene.gff3.sort.gff3.cds.fa"
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa"
     shell  : "python2.7 translation.py {input} {output}"

rule newgenes2uniprot:
     input  : "finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa.new_gene.fa"
     output : "blastout/newgene.na1.ev1e2.bo"
     params : uniprot=UNIPROT
     threads: 5 
     shell  : "blastp -query {input} -db {params.uniprot} -out {output} -outfmt 7 -num_threads {threads} -num_alignments 1 -evalue 1e-2"  

rule uniprot2annot:
     input  : "blastout/newgene.na1.ev1e2.bo"
     output : "finalout/newgene.annot"
     params : uniprot=UNIPROT
     shell  : "python2.7 uniprot_annot.py {input} {params.uniprot}; mv blastout/newgene.na1.ev1e2.bo.annot {output} "
     
