#!/bin/env python

fout = open('/g/data3/ba08/ctcf/ctctf_single_position_with_context_bed', 'w')
with open('/g/data3/ba08/ctcf/CTCF_motif_2000_regions.ed.bed') as fin:
    for line in fin:
        bits = line.strip().split('\t')
        chr=bits[0]
        orientation=bits[3]
        #if orientation == '+':
        position = -1001
        #elif orientation == '-':
        #    position = -1000
#        start_of_motif=int(bits[1]) - 1
#        end_of_motif=int(bits[2]) - 1
        start_of_motif=int(bits[1])
        end_of_motif=int(bits[2])
   
        motif=bits[4]
        if orientation == '+':
#            start_of_motif_inc = start_of_motif + 1
            start_of_motif_inc = start_of_motif - 1
            
            while start_of_motif_inc < end_of_motif:
                position += 1
                fout.write(str(chr) + '\t'+ str(start_of_motif_inc) + '\t' + str(start_of_motif) + '\t' + str(position) + '\t' + str(orientation) + '\t' + str(motif) + '\n')
                start_of_motif += 1
                start_of_motif_inc +=1
        elif orientation == '-':
            #end_of_motif_inc=int(bits[2]) - 1
            end_of_motif_inc=end_of_motif - 1
            while end_of_motif_inc > start_of_motif :
                position += 1
                fout.write(str(chr) + '\t'+ str(end_of_motif) + '\t' + str(end_of_motif_inc) + '\t' + str(position) + '\t' + str(orientation) + '\t' + str(motif) + '\n')
                end_of_motif -= 1
                end_of_motif_inc -=1
