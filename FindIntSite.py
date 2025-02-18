#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class FindIntSite:
    """
    A class to identify integration sites in sequencing data using BBMap and Clustal Omega.

    This class aligns sequencing reads to a reference sequence using BBMap, identifies soft-clipped reads,
    and uses Clustal Omega to align and analyze overhanging sequences. The results are saved in an output file.

    Attributes:
        input_1 (str): Path to the first input FASTQ file.
        input_2 (str): Path to the second input FASTQ file.
        inserted_seq (str): Path to the reference sequence file (inserted sequence).

    Methods:
        run_bbmap(): Aligns FASTQ files to the reference sequence using BBMap.
        parse_soft_clipped_reads(): Parses soft-clipped reads from the SAM file.
        create_clustal_input(f_list): Creates an input file for Clustal Omega.
        run_clustal_omega(): Runs Clustal Omega to align sequences.
        process_clustal_output(): Processes the output from Clustal Omega.
        find_overhanging_sequences(final_seq_list): Identifies overhanging sequences from the Clustal Omega output.
        get_unique_longest_sequences(seq_list): Filters unique and longest overhanging sequences.
        count_subsequence_overlaps(unique_longest_seqs, subsequences): Counts overlapping subsequences.
        generate_output(result_dict_5_end, result_dict_3_end): Generates the final output file.
        cleanup(): Removes temporary files.
        run(): Executes the entire workflow.
    """
        
    def __init__(self, input_1, input_2, inserted_seq):
        """
        Initializes the FindIntSite class with input file paths.

        Args:
            input_1 (str): Path to the first input FASTQ file.
            input_2 (str): Path to the second input FASTQ file.
            inserted_seq (str): Path to the reference sequence file.
        """
        self.input_1 = input_1
        self.input_2 = input_2
        self.inserted_seq = inserted_seq

    def run_bbmap(self):
        """Run bbmap to align fastq files to the reference sequence."""
        try:
            if self.input_2 is None:
                script = f"bbmap.sh minratio=0.2 in1={self.input_1} ref={self.inserted_seq} out=mapped.sam"
            else:
                script = f"bbmap.sh minratio=0.2 in1={self.input_1} in2={self.input_2} ref={self.inserted_seq} out=mapped.sam"
            os.system(script)
            logging.info("BBMap alignment completed successfully.")
        except Exception as e:
            logging.error(f"Error running BBMap: {e}")
            raise

    def parse_soft_clipped_reads(self):
        """Parse soft-clipped reads from the SAM file."""
        lines = []
        with open('mapped.sam', 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    lines.append(line.split())
        return [ls for ls in lines if "S" in ls[5]]

    def create_clustal_input(self, f_list):
        """Create input file for Clustal Omega."""
        ref_seq = ''
        with open(self.inserted_seq, "r") as i_file:
            for line in i_file:
                if not line.startswith('>'):
                    ref_seq += line.strip()
        with open("soft_reads.fasta", 'w') as s_file:
            s_file.write(">ref\n" + ref_seq + "\n")
            for i, ls in enumerate(f_list):
                s_file.write(f">seq{i+1}\n" + ls[-4] + "\n")

    def run_clustal_omega(self):
        """Run Clustal Omega to align sequences."""
        try:
            os.system("clustalo -i soft_reads.fasta -o aln.fasta -v")
            logging.info("Clustal Omega alignment completed successfully.")
        except Exception as e:
            logging.error(f"Error running Clustal Omega: {e}")
            raise
    def process_clustal_output(self):
        """Process the output from Clustal Omega."""
        with open('aln.fasta', 'r') as f:
            seq_str = f.read()
            seq_list = seq_str.split('>')
        final_seq_list = [''.join([l for l in seq if l in ["-", 'A', 'C', 'G', 'T']]) for seq in seq_list if seq]
        return final_seq_list
    
    def find_overhanging_sequences(self, final_seq_list: list):
        """Find overhanging sequences from the Clustal Omega output."""
        ref_seq = final_seq_list[0]
        start_index_5 = len(ref_seq) - len(ref_seq.lstrip('-'))
        start_index_3 = len(ref_seq) - len(ref_seq.rstrip('-'))

        final_seq_list.sort()
        ref_seq_index = final_seq_list.index(ref_seq)
        
        #Separating 5' and 3' end sequences.
        list_of_5_seqs = final_seq_list[ref_seq_index+1:]
        list_of_3_seqs = final_seq_list[:ref_seq_index]
        
        #Removing parts of the sequences that belong to the inserted sequence and cleaning the reminant seqs.
        list_of_5_seqs_cln = [seq[:start_index_5].replace("-", "") for seq in list_of_5_seqs]
        list_of_3_seqs_cln = [seq[-start_index_3:].replace("-", "") for seq in list_of_3_seqs]
        
        #Removing duplicate entries and sorting lists for further processing.
        list_of_5_seqs_cln = list(set(list_of_5_seqs_cln))
        list_of_5_seqs_cln.sort(key=lambda x: x[::-1], reverse=True)
        list_of_3_seqs_cln = list(set(list_of_3_seqs_cln))
        list_of_3_seqs_cln.sort(reverse=True)
        return list_of_5_seqs_cln, list_of_3_seqs_cln

    def get_unique_longest_sequences(self, seq_list: list):
        """Get unique longest overhanging sequences."""
        results = []
        for item in (seq_list):
            # Check if item is not already in results and not a substring of any sequence in results.
            if item not in results and all(item not in seq for seq in results):
                results.append(item)
        return results
    
    def count_subsequence_overlaps(self, unique_longest_seqs: list, subsequences: list):
        """
        This methods returns the metrics to to measure the confidence in the final returned results. It counts
        how many subsequences from unique_longest_seqs appear in each sequence in unique_longest_5_end as values
        to each unique sequence key in the result dictionary returned by the method. The higher 
        the number the more confidence in the results.

        Args:
            results_5_end (list): List of target sequences.
            list_of_5_seqs_cln (list): List of subsequences to search for.

        Returns:
            dict: A dictionary where keys are sequences from unique_longest_seqs and values are counts of matching 
            subsequences.
        """
        #Removing sequences less than 5 nucleotides to avoid aligning to more than one unique_longest_seq and giving false count.
        subsequences = list(filter(lambda item: len(item) >= 5,subsequences))
        overlap_dict = {}

        for target in unique_longest_seqs:
            count = 0
            # Count how many subsequences appear in the target sequence
            for subseq in subsequences:
                if subseq in target:
                    count += 1
            overlap_dict[target] = count

        return overlap_dict
    

    def generate_output(self, result_dict_5_end: dict, result_dict_3_end: dict):
        """Generate the final output file."""
        file_str = "Overhanging sequences at 5' end\n"
        for key, value in result_dict_5_end.items():
            file_str +=f"{key}: {value}\n"
        file_str += "\nOverhanging sequences at 3' end\n"
        for key, value in result_dict_3_end.items():
            file_str +=f"{key}: {value}\n"
        with open("BLAT_ready_seqs.txt", 'w') as s_file:
            s_file.write(file_str)
        return file_str

    def cleanup(self):
        """Remove temporary files."""
        os.remove('aln.fasta')
        os.remove('soft_reads.fasta')
        os.remove('mapped.sam')

    def run(self):
        """Run the entire workflow."""
        self.run_bbmap()
        f_list = self.parse_soft_clipped_reads()
        self.create_clustal_input(f_list)
        self.run_clustal_omega()
        final_seq_list = self.process_clustal_output()
        
        list_of_5_seqs_cln, list_of_3_seqs_cln = self.find_overhanging_sequences(final_seq_list)
        
        unique_longest_5_end = self.get_unique_longest_sequences(list_of_5_seqs_cln)
        unique_longest_3_end = self.get_unique_longest_sequences(list_of_3_seqs_cln)
        
        result_dict_5_end = self.count_subsequence_overlaps(unique_longest_5_end, list_of_5_seqs_cln)
        result_dict_3_end = self.count_subsequence_overlaps(unique_longest_3_end, list_of_3_seqs_cln) 
        
        output = self.generate_output(result_dict_5_end, result_dict_3_end)
        self.cleanup()
        return output
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find Integration site of inserted DNA sequences")
    parser.add_argument('-in1', '--input_1', type=str, required=True, help='The NGS fastq file with the inserted gene as file.fq')
    parser.add_argument('-in2', '--input_2', type=str, help='The NGS fastq file with the inserted gene as file.fq')
    parser.add_argument('-ref', '--inserted_seq', type=str, required=True, help='The inserted sequence in fasta format as file.fa')
    args = parser.parse_args()

    find_int_site = FindIntSite(args.input_1, args.input_2, args.inserted_seq)
    print("BLAT search ready sequences:\n", find_int_site.run())