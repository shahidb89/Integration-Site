#!/usr/bin/env python
# coding: utf-8

# In[42]:


import os
import argparse



parser = argparse.ArgumentParser(description = "Find Integration site of inserted DNA sequences")
parser.add_argument('-in1','--input_1', type=str, help='The NGS fastq file with the inserted gene as file.fq')
parser.add_argument('-in2','--input_2', type=str, help='The NGS fastq file with the inserted gene as file.fq')
parser.add_argument('-ref','--inserted_seq', type=str, help='The inserted sequence in fasta format as file.fa')
args = parser.parse_args()



# Run bbmap using command line script
def get_BLAT_seqs(input_1, input_2, inserted_seq):
    


    if input_2 == None:
        script = "bbmap.sh minratio=0.2 in1={} ref={} out=mapped.sam ".format(input_1, inserted_seq)
    else:
        script = "bbmap.sh minratio=0.2 in1={} in2={} ref={} out=mapped.sam ".format(input_1, input_2, inserted_seq)
    


    os.system(script)

    # Parse reads with "S" CIGAR string
    lines = []

    with open('mapped.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                pass
            else:
                lines.append(line.split())


    f_list = []

    for ls in lines:
        if "S" in ls[5]:
            f_list.append(ls)


    seq = ''

    with open('LTR.fa', "r") as LTR:
            for line in LTR:
                if not line.startswith('>'):
                    seq += line.strip()
    c = 1               
    with open("soft_reads.fasta", 'a') as s_file:
        s_file.write(">ref\n" + seq + "\n")
        for ls in f_list:
            s_file.write(">seq" + str(c) + "\n" + ls[-4] + "\n")
            c += 1


    os.system("clustalo -i soft_reads.fasta -o aln.fasta -v  ")    

    with open('aln.fasta', 'r') as f:
        seq_str = f.read()

    seq_list = seq_str.split('>')

    final_seq_list = []

    for seq in seq_list:
        temp = ''
        for l in seq:
            if l in "-" or l in 'A' or l in 'C' or l in 'G' or l in 'T' :
                temp +=l
        final_seq_list.append(temp)


    del final_seq_list[0]
    final_seq_list

    final_seq_list_sorted = sorted(final_seq_list)

    final_seqs_cln = []
    temp = ''

    for l in final_seq_list[0]:
        if l != '-':
            temp += l
    final_seqs_cln.append(temp)
    temp = ''

    for l in final_seq_list_sorted[0]:
        if l != '-':
            temp += l
    final_seqs_cln.append(temp)
    temp = ''
    for l in final_seq_list_sorted[-1]:
        if l != '-':
            temp += l
    final_seqs_cln.append(temp)



    # Capturing the length of overhanging sequences
    ltr_5 = 0
    ltr_3 = 0
    a = -1
    for n in final_seq_list[0]:
        if n in '-':
            ltr_5 += 1
        else:
            break
    c = -1       
    for n in range (len(final_seq_list[0])):
        if final_seq_list[0][c] in '-':
            ltr_3 += 1
            c -= 1
        else:
            break
    # Creating BLAT serach ready final list of sequences

    BLAT_ready_list = []

    BLAT_ready_list.append(final_seqs_cln[2][:ltr_5])
    BLAT_ready_list.append(final_seqs_cln[1][(len(final_seqs_cln[1])-ltr_3):])
    
    
    with open("BLAT_ready_seqs.txt", 'a') as s_file:
         s_file.write(">5'seq " +  BLAT_ready_list[0] + "\n" + ">3' seq " + BLAT_ready_list[1])
         
    
    os.remove('aln.fasta') 
    os.remove('soft_reads.fasta')
    os.remove('mapped.sam')
    
    return BLAT_ready_list


if __name__ == '__main__':
    print ("BLAT search ready 5' and 3' sequences consequently:", get_BLAT_seqs(args.input_1, args.input_2, args.inserted_seq))