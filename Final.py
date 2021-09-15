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


    os.system("clustalo -i soft_reads.fasta -o aln.fasta -v  ")    

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

    del final_seq_list[0]
    final_seq_list[0]

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

    final_seq_list.sort()
    for i, seq in enumerate(final_seq_list):
        if ref_seq in seq:
            ref_seq_index = i

    list_of_5_seqs = final_seq_list[ref_seq_index+1:]
    list_of_3_seqs = final_seq_list[:ref_seq_index]


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
        if seq_1 == seq_2:
            pass
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
    # Creating BLAT serach ready final list of sequences

    BLAT_ready_list = []
    BLAT_ready_list.append(results_5_end)
    BLAT_ready_list.append(results_3_end)
    
    file_str = ">5' overhanging seqs\n "
    for element in results_5_end:
        file_str += '>' + element + '\n'
        
    file_str += "\n>3' overhanging seqs\n "
    for element in results_3_end:
        file_str += '>' + element + '\n'
    
    
    with open("BLAT_ready_seqs.txt", 'w') as s_file:
         s_file.write(file_str)
         
    
    os.remove('aln.fasta') 
    os.remove('soft_reads.fasta')
    os.remove('mapped.sam')
    
    return BLAT_ready_list


if __name__ == '__main__':
    print ("BLAT search ready 5' and 3' sequences consequently:", get_BLAT_seqs(args.input_1, args.input_2, args.inserted_seq))