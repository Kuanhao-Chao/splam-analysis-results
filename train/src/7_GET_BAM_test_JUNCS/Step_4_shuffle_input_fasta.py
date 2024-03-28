# import random

# def shuffle_fasta(fasta_file, output_file):
#     """
#     Shuffle the order of sequences in a FASTA file.

#     Parameters:
#     fasta_file (str): Path to the input FASTA file.
#     output_file (str): Path to the output FASTA file with shuffled sequences.
#     """
#     # Read sequences and headers from the FASTA file
#     sequences = []
#     current_seq = ''
#     with open(fasta_file, 'r') as f:
#         for line in f:
#             if line.startswith('>'):
#                 if current_seq:
#                     sequences.append((header, current_seq))
#                     current_seq = ''
#                 header = line.strip()
#             else:
#                 current_seq += line.strip()
#         # Don't forget to save the last sequence
#         if current_seq:
#             sequences.append((header, current_seq))
    
#     # Shuffle the sequences
#     random.shuffle(sequences)
    
#     # Write the shuffled sequences to the output file
#     with open(output_file, 'w') as out:
#         for header, seq in sequences:
#             out.write(f"{header}\n")
#             out.write(f"{seq}\n")
import random
from Bio import SeqIO

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
SEQ_LEN = "800"
HALF_SEQ_LEN = int(SEQ_LEN) // 2
QUATER_SEQ_LEN = int(SEQ_LEN) // 4

# Example usage
output_dir = f'{project_root}train/results/BAM_TEST_junctions/INPUTS/'
input_fasta = f'{output_dir}{SEQ_LEN}bp/input_test.fa'
output_fasta = f'{output_dir}{SEQ_LEN}bp/input_test_shuffled.fa'

# Read the FASTA file
sequences = list(SeqIO.parse(input_fasta, "fasta"))

# Shuffle the sequences
random.shuffle(sequences)

# Write the shuffled sequences to a new file
with open(output_fasta, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta-2line")

