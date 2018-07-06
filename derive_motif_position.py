#!/bin/env python

fout = open('/g/data3/ba08/ctcf/ctctf_single_position_bed', 'w')
with open('/g/data3/ba08/ctcf/CTCF_motif.ed.bed') as fin:
    for line in fin:
        position = 0
        bits = line.strip().split('\t')
        chr=bits[0]
        orientation=bits[3]
        start_of_motif=int(bits[1]) - 1
        end_of_motif=int(bits[2]) - 1
        motif=bits[4]
        if orientation == '+':
            start_of_motif_inc = start_of_motif + 1
            while start_of_motif < end_of_motif:
                position += 1
                fout.write(str(chr) + '\t'+ str(start_of_motif) + '\t' + str(start_of_motif_inc) + '\t' + str(position) + '\t' + str(orientation) + '\t' + str(motif) + '\n')
                start_of_motif += 1
                start_of_motif_inc +=1
        elif orientation == '-':
            #end_of_motif_inc=int(bits[2]) - 1
            end_of_motif_inc=end_of_motif - 1
            while end_of_motif > start_of_motif :
                position += 1
                fout.write(str(chr) + '\t'+ str(end_of_motif_inc) + '\t' + str(end_of_motif) + '\t' + str(position) + '\t' + str(orientation) + '\t' + str(motif) + '\n')
                end_of_motif -= 1
                end_of_motif_inc -=1
