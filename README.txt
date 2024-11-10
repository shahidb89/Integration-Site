Here is a `README.md` file for the `FindIntSite.py` script, designed to provide a clear introduction, installation guidance, usage instructions, and output examples.

---

# **FindIntSite.py**

**FindIntSite.py** is a fast, automated Python tool for detecting integration sites of inserted DNA sequences in gene-editing and genetic engineering contexts. Leveraging the alignment power of BBMap and Clustal Omega, this tool processes whole genome sequencing (WGS) data, finds unique overhanging sequences (OVS), and prepares them for further search via BLAT to pinpoint integration sites at single-nucleotide precision.

### **Author**
- **Shahid HADI** (c) 2021
- **Contact**: [shahid_b89@hotmail.com](mailto:shahid_b89@hotmail.com)

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




