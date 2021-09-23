#!/usr/bin/env python
# coding: utf-8

# FindInSite is an automated python workflow to detect integration sites of inserted genes in gene edting procedures using whole genome sequencing approach. 

# Author: Shahid HADI
#Email: shahid_b89@hotmail.com



import os
import argparse



parser = argparse.ArgumentParser(description = "Find Integration site of inserted DNA sequences")
parser.add_argument('-in1','--input_1', type=str, help='The NGS fastq file with the inserted gene as file.fq')
parser.add_argument('-in2','--input_2', type=str, help='The NGS fastq file with the inserted gene as file.fq')
parser.add_argument('-ref','--inserted_seq', type=str, help='The inserted sequence in fasta format as file.fa')
args = parser.parse_args()



def get_BLAT_seqs(input_1, input_2, inserted_seq):
    """This function perfoms multiple tasks as following:
    1- It aligns fastq files obtained from Next Generation Sequencing NGS output to a referece sequence, 
     that is the    sequence inserted into the genome being  passed through NGS, using short aligner 
     program bbmap.
    2- Then, it selects all the soft-clipped aligend reads and puts them into a file with the reference 
     inserted sequence.
    3- The file is then passed in Clustal Omega so that all the soft_clipped sequence are aligned to the 
     reference sequence.
    4- The results fasta file of the Clustal Omega aligning is then processed to return the uniqe longest
     overhanging sequences and a number with each returned sequnece which indicates the number of other
     sequences that are both aligned to the referenece sequence and at the same time contain subsequences
     of thr returned sequence. The greater the number, the more confident the result, as longest 
     overhanging sequences have only a short part of them being aligned to the reference sequence.
     All the sequences are printed to the screen and a file named "BLAT_ready_seqs.txt".
     These sequences can be used in BLAT search to indictate the site of integration of the inserted 
     sequence.
    
    The function takes three inpute files. Two fastq files with .fq suffix(file.fq) and one fasta file 
    of the reference seq with .fa suffix (ref.fa)
    """

# Run bbmap using command line script
    if input_2 == None:
        script = "bbmap.sh minratio=0.2 in1={} ref={} out=mapped.sam ".format(input_1, inserted_seq)
    else:
        script = "bbmap.sh minratio=0.2 in1={} in2={} ref={} out=mapped.sam ".format(input_1, input_2, inserted_seq)
    


    os.system(script)

# Parse reads with "S" CIGAR string, thta is selecting soft-clipped reads from the sam file
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

# Creating the input file to Clustal Omega 
    ref_seq = ''

    with open('{}'.format(inserted_seq), "r") as i_file:
            for line in i_file:
                if not line.startswith('>'):
                    ref_seq += line.strip()
    c = 1               
    with open("soft_reads.fasta", 'w') as s_file:
        s_file.write(">ref\n" + ref_seq + "\n")
        for ls in f_list:
            s_file.write(">seq" + str(c) + "\n" + ls[-4] + "\n")
            c += 1

# Running the Clustal Omega program
    os.system("clustalo -i soft_reads.fasta -o aln.fasta -v  ")   

# Importing the results file from Clustal Omega alignment
    with open('aln.fasta', 'r') as f:
        seq_str = f.read()

        seq_list = seq_str.split('>')

    final_seq_list = []

    for seq in seq_list:
        temp = ''
        for l in seq:
            if l in ["-", 'A', 'C', 'G', 'T'] :
                temp +=l
        final_seq_list.append(temp)
# Cleaning the list
    del final_seq_list[0]

# In the fasta file returned by ClustalO, all the sequences have the same length of characteres with empty spaces being given a (-), final_seq_list[0] is the reference sequence, here we are capturing the index were the soft-clipped sequences divide into to distinct subsequences, one aligned to the reference genome, and another aligned to the host genome. We do this for both ends of the reference sequence. Here the two ends of the reference sequence behave differently! So, each should be handled carefully! 
    start_index_5 = 0
    start_index_3 = 0

    for i in final_seq_list[0]:
        if i in "-":
            start_index_5 +=1
        else:
            break

    c = -1       
    for n in range (len(final_seq_list[0])):
        if final_seq_list[0][c] in '-':
            start_index_3 += 1
            c -= 1
        else:
            break
# The final_seq_list is sorted so that the reference sequence comes in the middle and divides soft-clipped sequences to two. The first part are the sequences that align to the 3' (righ side) end of the reference sequences(list_of_3_seqs), and the last part are the sequences that align to the 5' (left side) of the reference sequence (list_of_5_seqs).
    final_seq_list.sort()
    for i, seq in enumerate(final_seq_list):
        if ref_seq in seq:
            ref_seq_index = i

    list_of_5_seqs = final_seq_list[ref_seq_index+1:]
    list_of_3_seqs = final_seq_list[:ref_seq_index]

# Cleaning both lists
    list_of_5_seqs_cln = []
    list_of_3_seqs_cln = []

    for i in list_of_5_seqs:
        i = i[:start_index_5].replace("-","")
        list_of_5_seqs_cln.append(i)

    def sort_5_seqs(seq):
        return seq[::-1]

    list_of_5_seqs_cln = list(set(list_of_5_seqs_cln))
    list_of_5_seqs_cln_sorted = sorted(list_of_5_seqs_cln, key=sort_5_seqs)

    list_of_5_seqs_cln_sorted.append("")

    for i in list_of_3_seqs:
        i = i[-start_index_3:].replace("-",'')
        list_of_3_seqs_cln.append(i)

    list_of_3_seqs_cln = list(set(list_of_3_seqs_cln))
    list_of_3_seqs_cln_sorted = sorted(list_of_3_seqs_cln) 

    list_of_3_seqs_cln_sorted.append("")

    results_3_end = []
# The following code returns a list of unique longet overhanging sequences
    a = 0
    b = 1
    c = 2


    while a < (len(list_of_3_seqs_cln_sorted)-2):
        seq_1 = list_of_3_seqs_cln_sorted[a]
        seq_2 = list_of_3_seqs_cln_sorted[b]
        seq_3 = list_of_3_seqs_cln_sorted[c]
        if seq_1 in seq_2 and seq_2 not in seq_1 and len(seq_1)<len(seq_2)>len(seq_3)and seq_2.startswith(seq_1):
            results_3_end.append(seq_2)
            a+=1
            b+=1
            c+=1
        else:
            pass

            a+=1
            b+=1
            c+=1
            
    results_5_end = []
    a = 0
    b = 1
    c = 2
    while a < (len(list_of_5_seqs_cln_sorted)-2):
        seq_1 = list_of_5_seqs_cln_sorted[a]
        seq_2 = list_of_5_seqs_cln_sorted[b]
        seq_3 = list_of_5_seqs_cln_sorted[c]
        if (seq_1 in seq_2) and (seq_2 not in seq_1) and len(seq_1)<len(seq_2)>len(seq_3) and seq_2.endswith(seq_1):
            results_5_end.append(seq_2)
            a+=1
            b+=1
            c+=1
        else:
            pass
            a+=1
            b+=1
            c+=1
# The following code captures the number of intermediatry sequences between a unique longest sequence returned and the refereces sequence. 

    l5_subseq_count = []
    l5 = list_of_5_seqs_cln_sorted
    for i in  results_5_end:
        a = 0
        for n in l5:
            if i != n:
                a += 1
            else:
                l5_subseq_count.append(a)
                l5 = l5[a+1:]
                
    l3_subseq_count = []
    l3 = list_of_3_seqs_cln_sorted
    for i in  results_3_end:
        a = 0
        for n in l3:
            if i != n:
                a += 1
            else:
                l3_subseq_count.append(a)
                l3 = l3[a+1:]
# Building the returned string and file
    file_str = "5' overhanging seqs\n"
    for (element,number) in zip(results_5_end, l5_subseq_count):
        file_str += '>' + element + '\t' + str(number) + '\n'
        
    file_str += "\n3' overhanging seqs\n"
    for (element, number) in zip(results_3_end, l3_subseq_count):
        file_str += '>' + element + '\t' + str(number) + '\n'
    
    
    with open("BLAT_ready_seqs.txt", 'w') as s_file:
         s_file.write(file_str)
         


    os.remove('aln.fasta') 
    os.remove('soft_reads.fasta')
    os.remove('mapped.sam')
    
    return file_str


if __name__ == '__main__':
    print ("BLAT search ready 5' and 3' sequences consequently:\n", get_BLAT_seqs(args.input_1, args.input_2, args.inserted_seq))