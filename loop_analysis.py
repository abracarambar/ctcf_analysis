#!/bin/env python
import subprocess
import os
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '/g/data2/gd7/software/maegau/miniconda2/envs/localpython')

#original bash commands to work out the steps to take
#grep -v NA Looplist_with_motifs_ed.txt | awk 'BEGIN{OFS="\t"; col13=0} {if ($9 == "p" && $12 == "n") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"3_CTCF_motif",$1,$10,$11,"loop_"col13,$1,$8,$10; col13+=1}' > Looplist_with_3motifs_txt
#grep -v NA Looplist_with_motifs_ed.txt | awk 'BEGIN{OFS="\t"; col13=0} {if ($9 == "p" && $12 == "n") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"5_CTCF_motif",$1,$7,$8,"loop_"col13,$1,$8,$10; col13+=1}' > Looplist_with_5motifs_txt
#bedtools intersect -wa -wb -b  ucsc.bed -a Looplist_txt | cut -f1-4,10 | sort | uniq > Genes_within_loops_txt
#grep -v NA Looplist_with_motifs_ed.txt | awk 'BEGIN{OFS="\t"; col13=0} {if ($9 == "p" && $12 == "n") print$1,$8,$10,"loop_"col13; col13+=1}' | sort -k1,1 -k2,2 -V > Looplist_txt


fout = open('/g/data3/ba08/ctcf/motifs_with_loops.bed', 'w')
fout2 = open('/g/data3/ba08/ctcf/loops.bed', 'w')
with open('/g/data3/ba08/ctcf/Looplist_with_motifs_ed.txt') as fin:
    loop_count = 0
    for line in fin:
        if "NA" not in line and "chr1" not in line:
            bits = line.strip().split('\t')
            chr=bits[0]
            start_of_5_ctcf_motif=int(bits[6])
            end_of_5_ctcf_motif=int(bits[7])
            orientation_of_5_ctcf_motif=str(bits[8])
            orientation_of_3_ctcf_motif=str(bits[11])
            start_of_3_ctcf_motif=int(bits[9])
            end_of_3_ctcf_motif=int(bits[10])
            start_of_loop=end_of_5_ctcf_motif
            end_of_loop=start_of_3_ctcf_motif
            loop_count += 1
            if orientation_of_5_ctcf_motif == "p" and orientation_of_3_ctcf_motif == "n" and start_of_5_ctcf_motif != "NA" and start_of_3_ctcf_motif != "NA":
		fout.write(str(chr) + '\t' + str(start_of_5_ctcf_motif) + '\t' + str(end_of_5_ctcf_motif) + '\t+\t5_CTCF\tloop_' + str(loop_count) + '\n')
	    	fout.write(str(chr) + '\t' + str(start_of_3_ctcf_motif) + '\t' + str(end_of_3_ctcf_motif) + '\t+\t3_CTCF\tloop_' + str(loop_count) + '\n')
                fout2.write(str(chr) + '\t' + str(start_of_loop) + '\t' + str(end_of_loop) + '\t+\tloop_' + str(loop_count) + '\n')

fout2.close()
fout.close()

#sort loop bed file
loops_bed='/g/data3/ba08/ctcf/loops.bed'
sorted_loops_bed = '/g/data3/ba08/ctcf/loops_sorted.bed'
cmd1 = ['sort', '-k1,1', '-k2,2n' , '-V', loops_bed, '-o', sorted_loops_bed]
stdout = open('/g/data3/ba08/ctcf/error_log',"wb")
subprocess.check_call(cmd1, stdout=stdout)
os.remove(loops_bed)

motifs_bed= '/g/data3/ba08/ctcf/motifs_with_loops.bed'
sorted_motifs_bed = '/g/data3/ba08/ctcf/motifs_with_loops_sorted.bed'
cmd2 = ['sort', '-k1,1', '-k2,2n', '-V', motifs_bed, '-o', sorted_motifs_bed]
subprocess.check_call(cmd2, stdout=stdout)
os.remove(motifs_bed)

#extract genes present in each loop
#output bed file with coordinate of loops and gene within loop on individual rows

ucsc_bed = '/g/data3/ba08/ctcf/ucsc.bed'
bedtools_cmd = ['bedtools', 'intersect', '-wa', '-wb', '-b', ucsc_bed, '-a', sorted_loops_bed, '|', 'cut', '-f1-5,11', '|', 'sort', '|', 'uniq', '>', 'Genes_within_loops_txt']
bedtools_cmd_str = " ".join(bedtools_cmd)
bedtools_cmd_proc = subprocess.Popen(bedtools_cmd_str, shell=True)
bedtools_cmd_proc.communicate()
sample_list = [183410, 184577, 193958, 200971, 321773, 33432, 34366, 34934, 35562, 35649, 35818, 38532, 4699, 48585, 9120]
#sample_list = [183410, 184577]
full_gene_list=[]
full_gene_count=0
full_motif_count=0
full_motif_list=[]
new_dict={}
summary_out=open('/g/data3/ba08/ctcf/summary.txt', 'w')
summary_out.write('Sample' + '\t' + 'Mutated_Motifs_Count' + '\t' + 'Associated_Loop_Count' + '\t' + 'Within_Loop_Gene_Count' + '\n')
for sample in sample_list:
    vcf = '/g/data2/gd7/software/maegau/vcf2maf/tests/' + str(sample) + '.strelka.filtered.pass.vcf'
    bedtools_cmd2 = ['bedtools', 'intersect', '-wa', '-wb', '-b', vcf, '-a', sorted_motifs_bed, '>', '/g/data3/ba08/ctcf/overlap_between_' + str(sample) + '_and_motifs.txt']
    bedtools_cmd2_str = " ".join(bedtools_cmd2)
    bedtools_cmd2_proc = subprocess.Popen(bedtools_cmd2_str, shell=True)
    bedtools_cmd2_proc.communicate()
    loop_list=[]
    loop_count=0
    motif_count=0


    motif_list=[]
    with open('/g/data3/ba08/ctcf/overlap_between_' + str(sample) + '_and_motifs.txt') as fin:
        for line in fin:
            bits = line.strip().split('\t')
            loops_affected = bits[5]
            motif_affected= bits[4]
            motif_loop_affected ='%s-%s' % (loops_affected, motif_affected)
            #print motif_loop_affected
            if motif_loop_affected not in motif_list:
                motif_list.append(motif_loop_affected)
                motif_count += 1
            if motif_loop_affected not in full_motif_list:
                full_motif_list.append(motif_loop_affected)
                full_motif_count += 1
            if loops_affected not in loop_list:
                loop_list.append(loops_affected)
                loop_count += 1
    print 'motif count for', sample, motif_count
    print 'loop count for', sample, loop_count
    #print loop_list
    gene_count=0
    gene_list=[]
    for affected_loop in loop_list:
        with open('/g/data3/ba08/ctcf/Genes_within_loops_txt') as fin:
            for line in fin:
                bits = line.strip().split('\t')
                loop=bits[4]
                gene=bits[5]
                if affected_loop == loop:
                    #print gene
                    if gene not in full_gene_list:
                        full_gene_list.append(gene)
                        full_gene_count += 1
                    if gene not in gene_list:
                        gene_list.append(gene)
                        gene_count += 1
    new_dict[sample] = gene_list
    print 'gene count for', sample, gene_count
    #print gene_list 
    summary_out.write(str(sample) + '\t' + str(motif_count) + '\t' + str(loop_count) + '\t' + str(gene_count) + '\n')
print 'total gene count:', full_gene_count
print 'full motif count:', full_motif_count
#print new_dict

df = pd.DataFrame(full_gene_list, columns=['Genes'])
for sample in sample_list:
    print 'Adding genes occurring within loops with mutated CTCFs in their anchor for sample', sample
    gene_list_for_sample=new_dict[sample]
    #print gene_list_for_sample
    column=[]
    for gene in full_gene_list:
        if gene in gene_list_for_sample:
            column.append('1')
        else:
            column.append('0')
    se = pd.Series(column)
    df[sample] = se.values
#df = df.set_index('Genes')
#df.to_csv('/g/data3/ba08/ctcf/gene_table.txt', sep='\t')

csvfile = pd.read_csv('pancancer_gene_list.txt', header=0)
#print csvfile
#csvfile.index.name = 'Genes'
print csvfile
merged_df = pd.merge(df, csvfile, on='Genes')
merged_df.to_csv('/g/data3/ba08/ctcf/gene_table_nano.txt', sep='\t')
print merged_df
df = df.set_index('Genes')
df.to_csv('/g/data3/ba08/ctcf/gene_table.txt', sep='\t')
