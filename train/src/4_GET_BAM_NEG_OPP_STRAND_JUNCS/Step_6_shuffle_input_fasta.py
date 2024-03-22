import random

def shuffle_fasta(fasta_file, output_file):
    """
    Shuffle the order of sequences in a FASTA file.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    output_file (str): Path to the output FASTA file with shuffled sequences.
    """
    # Read sequences and headers from the FASTA file
    sequences = []
    current_seq = ''
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences.append((header, current_seq))
                    current_seq = ''
                header = line.strip()
            else:
                current_seq += line.strip()
        # Don't forget to save the last sequence
        if current_seq:
            sequences.append((header, current_seq))
    
    # Shuffle the sequences
    random.shuffle(sequences)
    
    # Write the shuffled sequences to the output file
    with open(output_file, 'w') as out:
        for header, seq in sequences:
            out.write(f"{header}\n")
            out.write(f"{seq}\n")
    

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 

SEQ_LEN = "800"
HALF_SEQ_LEN = int(SEQ_LEN) // 2
QUATER_SEQ_LEN = int(SEQ_LEN) // 4

# Example usage
for threshold in range(2, 101):
    output_dir = f'{project_root}train/results/Negs/Neg_{threshold}/INPUTS/'
    input_fasta = f'{output_dir}{SEQ_LEN}bp/input_neg_{threshold}.fa'
    output_fast = f'{output_dir}{SEQ_LEN}bp/input_neg_{threshold}_shuffled.fa'
    shuffle_fasta(input_fasta, output_fast)
