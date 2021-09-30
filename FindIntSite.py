#!/usr/bin/env python
# coding: utf-8

# FindIntSite is an automated python workflow to detect integration sites of inserted genes in gene edting procedures using whole genome sequencing approach. 

#Copyright: (c) Shahid HADI, 2021 Author: Shahid HADI Contact: shahid_b89@hotmail.com 



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
     overhanging sequences and teo numbers with each returned sequnece. The first number indicates the number 
     of ovelapping sequences between the returned sequence and the inserted sequence, while the second number         
     represents the number of the nucelotides in the returned sequence that align to the reference sequence.
     If the second number, the number of nucleotides in each sequence that map to the reference sequence, is small, 
     i.e. < 20, and the first number is zero or close to zero, then the confidence in the returned sequence is low and 
     needs further inspection.
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
# Capturing the length of the mapped reads that will be used to calculate the number of nucleotides of each returned 
# sequence that map to the reference sequence
    len_of_mapped_reads = len(f_list[1][-4])

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

# In the fasta file returned by ClustalO, all the sequences have the same length of characteres with empty 
# spaces being given a (-), final_seq_list[0] is the reference sequence, here we are capturing the index 
# were the soft-clipped sequences divide into to distinct subsequences, one aligned to the reference 
# genome, and another aligned to the host genome. We do this for both ends of the reference sequence. 
# Here the two ends of the reference sequence behave differently! So, each should be handled carefully! 
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
# The final_seq_list is sorted so that the reference sequence comes in the middle and divides 
# soft-clipped sequences to two. The first part are the sequences that align to the 3' (righ side) 
# end of the reference sequences(list_of_3_seqs), and the last part are the sequences that align to 
# the 5' (left side) of the reference sequence (list_of_5_seqs).
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
        elif seq_1 not in seq_2 and seq_2 not in seq_3:
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
        elif seq_1 not in seq_2 and seq_2 not in seq_3:
            results_5_end.append(seq_2)
            a+=1
            b+=1
            c+=1
        else:
            pass
            a+=1
            b+=1
            c+=1
# The following code captures the number of intermediatry sequences between a unique longest sequence 
# returned and the refereces sequence. 

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
                
# Capturing the number of nucleotides of each returned sequence that map to the reference sequence

    no_of_n_map_to_ref_5 = []
    no_of_n_map_to_ref_3 = []

    for m in results_5_end:
        y = len_of_mapped_reads - len(m)
        no_of_n_map_to_ref_5.append(y)

    for m in results_3_end:
        y = len_of_mapped_reads - len(m)
        no_of_n_map_to_ref_3.append(y)


# Building the returned string and file
    file_str = "Overhanging sequences at start\n"
    for (element,number, count) in zip(results_5_end, l5_subseq_count,  no_of_n_map_to_ref_5 ):
        file_str += '>' + element + '\t' + str(number) + '\t' + str(count) + '\n'
        
    file_str += "\nOverhanging sequences at end\n"
    for (element, number,count) in zip(results_3_end, l3_subseq_count, no_of_n_map_to_ref_3):
        file_str += '>' + element + '\t' + str(number) + '\t' + str(count) + '\n'
    
    
    with open("BLAT_ready_seqs.txt", 'w') as s_file:
         s_file.write(file_str)
         


    os.remove('aln.fasta') 
    os.remove('soft_reads.fasta')
    os.remove('mapped.sam')
    
    return file_str


if __name__ == '__main__':
    print ("BLAT search ready sequences:\n", get_BLAT_seqs(args.input_1, args.input_2, args.inserted_seq))