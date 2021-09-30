FindIntSite.py

Copyright: (c) Shahid HADI, 2021 Author: Shahid HADI Contact: shahid_b89@hotmail.com

FindIntSite.py is a handy, fast, easy to use, fully automated and all-purpose tool that can be used to detect any kind of integration sites in any context. The main concept is mapping the whole genome sequencing (WGS) data against the inserted sequence and then returning overhanging sequences (OVS). It works as an automated workflow that arranges a traffic between third party aligner programs, BBmap (Bushnell, 2015) and Clustal omega (Sievers et al., 2014), and our Python script Then, the overhanging sequences are searched for using BLAT to find the exact integration post at one nucleotide resolution.


FindIntSite.py takes takes a max to two WGS files in FASTQ format (WGS_file_1.fq, WGS_file_2.fq), and one file that contains the inserted sequence in FASTA format (insertted_seq.fa)

The command line input syntax is as follows:
python  FindIntSite.py -in1 WGS_file_1.fq -in2 WGS_file_2.fq -ref inserted_seq.fa 

Note that the second WGS input file is optional. So that the following command will work too:

python  FindIntSite.py -in1 WGS_file_1.fq -ref inserted_seq.fa 

A text file is returned that contains all the OVS being divided into two groups (Figure 4). OVS to the start of the reference (inserted) sequence, and OVS to the of the reference sequence. Each sequence is followed by two numbers. The first one represents the number of overlapping sequences between the returned sequence and the reference sequence and the second one represents the number of nucleotides in each sequence that map to the reference sequence.

Here is how  typical results file looks like with 8 insertion sites:

BLAT search ready sequences:git 
Overhanging sequences at start
>GTATGCTATCGGGCATGCGCGCTCATATTCTACGTAGGGATAATAATTATCAGTTTATTGTGACAA	17	34
>AACTTTCGTCCATTGCAGGTACCCCCATCGAGGCTGGGCCCTTCGATACAGTGGTGCTTTTTA	18	37
>TTCGGCCTGGAAAGCAAAAACTTTATGAAAAGGTGTAACCCGCCATCCATGGAAAAACCAGGCGC	18	35
>GCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGC	20	34
>TTGTCCGGCAATCCGCCAGTTGTGAATACCGCCCGCATGTCCGGTGCTTTTCAGCCCCAGTTTCCG	16	34
>ATATTTCATGTAAGGTTCATTATGAAATTAACTAAACTTGTACTTGAAAAT	8	49
>CTCATCAGGCTTTTTCACTGGTTTAGCGTTGATCAGACTGTTCATTCCCGCAGCAGTAGTT	18	39
>CAAAAATGGACAAGACTTGATCCTCTTCCATATTCACAGAGAATATTGAATATTGCACATGATTTT	18	34

Overhanging sequences at end
>AAATTAACTATAACTGATTTGCCATGATGCCCTTTGCAGCCGGAAAGGGCATCTCTTTCAGAAA	18	36
>CCTATAACGAAAGAAAGTGTTATTGTAGTTCATAATAGCATGAAGCTATACACTCAAATTCGTAA	15	35
>CGCCGAGTAGCTGGTCAGCAGCCACTGAGTCAGGATGCCATCTTTAATAATATCGCGACGCTCG	14	36
>CGGAGATTGGCGATGCGTTTGGTGGTCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATC	19	35
>GTTTGCGGAGTAATGTCTCGCTCAACGCGCGGTGCCGTTTCCTGTAATTCGTCAGGGGTGTA	21	38
>TGGTTATCAGATACTAAAAAGTTATTCATTATCGCCCGGGATGTCTTCTTGTCGGGAATTACCC	18	36
>TTAATAGTATAAATTACGGCTGCTATCTTGATAATACCGATATTAATTGGTGCAGCATTGTTTAT	19	35
>TTCAGATCTGTTAGAGAGCGGCAGGAAATTGATTTTGCACCTGTAACATTACTATTGGGG	20	40


