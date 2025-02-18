# **FindIntSite.py**
Updated Feb 2025

**FindIntSite.py** is a fast, automated Python tool for detecting integration sites of inserted DNA sequences in gene-editing and genetic engineering contexts. Leveraging the alignment power of BBMap and Clustal Omega, this tool processes whole genome sequencing (WGS) data, finds unique overhanging sequences (OVS), and prepares them for further search via BLAT to pinpoint integration sites at single-nucleotide precision.

### **Author**
- **Shahid HADI** 
- **Contact**: [shahid_b89@hotmail.com](mailto:shahid_b89@hotmail.com)

### Description

Understanding the integration sites of inserted genes is crucial in genetic engineering, especially in gene therapy applications where viral vectors are often used. These vectors can randomly integrate into transcriptionally active regions of the host genome, potentially disrupting important genes and causing unintended consequences. FindIntSite is designed to address this challenge by providing an automated and efficient method for detecting the precise locations of these insertions.

The process begins when a therapeutic gene, referred to as the "inserted gene," is introduced into a host genome with a defective gene. After the gene insertion, Whole Genome Sequencing (WGS) is performed on the modified genome. From this point, **FindIntSite** takes over to process the data.

Using **BBMap**, **FindIntSite** aligns the inserted sequence with the host genome's WGS data to identify the insertion points. This alignment generates **Overhanging Sequences (OVS)**, which contain portions of both the host genome and the inserted sequence at each integration site.

Next, **FindIntSite** separates the **OVS** into two distinct sequences: one representing the host genome and the other the inserted sequence. This separation is achieved by aligning the **OVS** with the inserted sequence using **Clustal Omega**.

**FindIntSite** then processes the host genome sequences to identify unique insertion sites, returning the longest sequences at each site and eliminating shorter, redundant sequences. The final output includes the identified integration sites, with details on whether the inserted sequence lies on the 5' or 3' end of the returned host genome sequence.

The identified integration sites are then ready for further analysis using **BLAT** to pinpoint the exact location of the insertion within the genome.


## **Features**
- **Automated Workflow**: Orchestrates alignment and sequence processing through BBMap and Clustal Omega.
- **Detects Unique Integration Sites**: Identifies and characterizes overhanging sequences around integration sites.
- **Detailed Output**: Provides sequences with nucleotide counts and overlap information, facilitating BLAT search for integration site analysis.

## **Dependencies**
To use **FindIntSite**, install the following dependencies:

- **BBMap**: [Download here](http://sourceforge.net/projects/bbmap/) (Requires Java).
- **Clustal Omega**: [Available on Anaconda](https://anaconda.org/bioconda/clustalo).

Alternatively you can download dependecies from the terminal by wrtiring the followning commands:
```bash
conda install bbmap
```
```bash
conda install clustalo
```

## **Installation**

1. Install **BBMap** and **Clustal Omega** as per their respective instructions.
2. Ensure Python 3 is installed.

## **Usage**

Run `FindIntSite.py` with the command line, specifying up to two WGS FastQ files and a single FASTA file containing the inserted sequence.

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

The script generates a text file, `BLAT_ready_seqs.txt`, with identified overhanging sequences at both the 5' and 3' ends of the inserted sequence. Each entry includes the sequence and a number separated by colon. The number is the number of the overlapping sequences between the give sequence and the inserted sequence. It can be used as a metric to measure the credibility of the results. The higher the number, the more reliable the reurned sequence is. This information can be directly used in BLAT for integration site localization.

### **Example Output**

```plaintext
BLAT search ready sequences:
 Overhanging sequences at 5' end
GTTTGCTCATCAGGCTTTTTCACTGGTTTAGCGTTGATCAGACTGTTCATTCCCGCAGCAGTAGTT: 16
TTGTCCGGCAATCCGCCAGTTGTGAATACCGCCCGCATGTCCGGTGCTTTTCAGCCCCAGTTTCCG: 17

Overhanging sequences at 3' end
GTTTGCGGAGTAATGTCTCGCTCAACGCGCGGTGCCGTTTCCTGTAATTCGTCAGGGGTGTAAACA: 22
CGCCGAGTAGCTGGTCAGCAGCCACTGAGTCAGGATGCCATCTTTAATAATATCGCGACGCTCGGT: 16
```

## **Explanation of Output Values**
Each sequence line in `BLAT_ready_seqs.txt` has:
1. **Overhanging Sequence**: Unique sequence at either end of the inserted DNA.
2. **Overlap Count**: Number of overlapping reads with the inserted sequence.


### **References**
- **BBMap**: Bushnell B. (2015). BBMap: a fast, accurate, splice-aware aligner. 
- **Clustal Omega**: Sievers, F. et al. (2014). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.

## **License**
This project is licensed under the MIT License.




