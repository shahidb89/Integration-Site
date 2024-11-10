# **FindIntSite.py**

**FindIntSite.py** is a fast, automated Python tool for detecting integration sites of inserted DNA sequences in gene-editing and genetic engineering contexts. Leveraging the alignment power of BBMap and Clustal Omega, this tool processes whole genome sequencing (WGS) data, finds unique overhanging sequences (OVS), and prepares them for further search via BLAT to pinpoint integration sites at single-nucleotide precision.

### **Author**
- **Shahid HADI** (c) 2021
- **Contact**: [shahid_b89@hotmail.com](mailto:shahid_b89@hotmail.com)

## **Description**
It is important to have tools for transgene insertion site discovery, as many gene therapy clinical trials use viral vectors that have the tendency to insert randomly into active transcriptionally active regions. This can be an issue if the function of important genes is disrupted. For a better comprehension of what is going on under the hood and to better understand the importance and scope of FindIntSite, we will go through an entire procedure that begins with gene insertion and ends with the discovery of insertion sites.
1. A therapeutic gene, we call it inserted gene, is inserted into a host genome that has a defective gene.
2. Whole Genome Sequencing, WGS, is performed for the genome after the insertion of the therapeutic (inserted) gene. From this step onwards, **FindIntSite** starts to do an amazing job.
3. Using **BBMap**, **FindIntSite** maps the inserted sequence onto the WGS results of the host genome to find where the inserted gene has been inserted across the host genome. This mapping returns multiple sequences that are known as Overhanging Sequences, OVS. These are sequences that contain parts of the host genome at each insertion site and parts of the inserted sequence.
4. **FindIntSite** further processes the OVSes by splitting them into two distinct sequences, separating the inserted genome part from the host genome part. This is done by mapping the OVSes onto the inserted sequence using **Clustal Omega**.
5. Then, **FindIntSite** further processes host genome sequences to identify unique insertion sites and return the longest sequence at each insertion site while getting rid of shorter sequences that belong to the same unique insertion site.
6. **FindIntSite** then returns as many insertion sites that have been captured by the WGS. It is important to mention that:
    - **FindIntSite** can return multiple insertion sites. As many as present in the WGS and captured by BBMap.
    - **FindIntSite** indicates where the inserted sequence lies relative to the returned sequence., i.e does the inserted sequence lie to the 5' end or 3' end of the returned host sequence?
7. Returned sequences are ready to be used for BLAT search.
## **Features**
- **Automated Workflow**: Orchestrates alignment and sequence processing through BBMap and Clustal Omega.
- **Detects Unique Integration Sites**: Identifies and characterizes overhanging sequences around integration sites.
- **Detailed Output**: Provides sequences with nucleotide counts and overlap information, facilitating BLAT search for integration site analysis.

## **Dependencies**
To use **FindIntSite**, install the following dependencies:

- **BBMap**: [Download here](http://sourceforge.net/projects/bbmap/) (Requires Java).
- **Clustal Omega**: [Available on Anaconda](https://anaconda.org/bioconda/clustalo).

## **Installation**

1. Install **BBMap** and **Clustal Omega** as per their respective instructions.
2. Ensure Python 3 is installed.

## **Usage**

Run `FindIntSite.py` with the command line, specifying up to two WGS FastQ files (for paired-end sequencing) and a single FASTA file containing the inserted sequence.

```bash
python FindIntSite.py -in1 WGS_file_1.fq -in2 WGS_file_2.fq -ref inserted_seq.fa
```

**Command-Line Arguments**:
- `-in1`, `--input_1`: First WGS FastQ file.
- `-in2`, `--input_2`: Optional second WGS FastQ file (for paired-end reads).
- `-ref`, `--inserted_seq`: FASTA file containing the inserted sequence.

If only one WGS file is used, the script works with a single FastQ input:

```bash
python FindIntSite.py -in1 WGS_file_1.fq -ref inserted_seq.fa
```

## **Output**

The script generates a text file, `BLAT_ready_seqs.txt`, with identified overhanging sequences at both the start and end of the inserted sequence. Each entry includes the sequence, the number of overlapping sequences with the reference, and the number of nucleotides mapping to the reference. This information can be directly used in BLAT for integration site localization.

### **Example Output**

```plaintext
BLAT search ready sequences:
Overhanging sequences at start
>GTATGCTATCGGGCATGCGCGCTCATATTCTACGTAGGGATAATAATTATCAGTTTATTGTGACAA	17	34
>AACTTTCGTCCATTGCAGGTACCCCCATCGAGGCTGGGCCCTTCGATACAGTGGTGCTTTTTA	18	37
...

Overhanging sequences at end
>AAATTAACTATAACTGATTTGCCATGATGCCCTTTGCAGCCGGAAAGGGCATCTCTTTCAGAAA	18	36
>CCTATAACGAAAGAAAGTGTTATTGTAGTTCATAATAGCATGAAGCTATACACTCAAATTCGTAA	15	35
...
```

## **Explanation of Output Values**
Each sequence line in `BLAT_ready_seqs.txt` has:
1. **Overhanging Sequence**: Unique sequence at either end of the inserted DNA.
2. **Overlap Count**: Number of overlapping reads with the reference sequence.
3. **Mapped Nucleotides**: Number of nucleotides mapped to the reference, indicating alignment accuracy.

### **References**
- **BBMap**: Bushnell B. (2015). BBMap: a fast, accurate, splice-aware aligner. 
- **Clustal Omega**: Sievers, F. et al. (2014). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.

## **License**
This script and documentation are copyrighted to Shahid HADI, 2021.




