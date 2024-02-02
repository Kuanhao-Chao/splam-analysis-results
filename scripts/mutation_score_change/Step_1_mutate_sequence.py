from Bio import SeqIO
import sys, os
import random

type = sys.argv[1]
target = "pos_MANE"
project_root = "../../train/src/INPUTS/800bp/"
output_dir = f"../../results/mutation/{target}/seq/"

def mutate_sequence(sequence, position, new_nucleotide):
    """Mutate the sequence at the specified position with the new nucleotide."""
    return sequence[:position] + new_nucleotide + sequence[position + 1:]

def random_mutation(original_nucleotide):
    nucleotides = ['A', 'T', 'C', 'G']
    nucleotides.remove(original_nucleotide)
    return random.choice(nucleotides)


def main(fasta_file):
    COUNTER = 5000
    counter = 0
    NO_MUTATION = True

    # If NO_MUTATION is True, the original sequence will be written to the output file
    if NO_MUTATION:
        no_mutated_filename = f"{output_dir}original.fasta"
        print("processing original sequence...")
        with open(no_mutated_filename, "w") as output_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                original_sequence = str(record.seq)
                sequence_id = record.id
                print(f"\tProcessing sequence {sequence_id}...")
                if "N" in original_sequence:
                    continue
                output_file.write(f">{sequence_id}\n")
                output_file.write(original_sequence + "\n")
                counter += 1
                if counter >= COUNTER:
                    break
        return 0
    
    # If NO_MUTATION is False, the original sequence will be mutated and written to the output file
    nucleotides = ['A', 'T', 'C', 'G']
    for position in range(800):
        mutated_filename = f"{output_dir}pos_{position}.fasta"
        with open(mutated_filename, "w") as output_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                original_sequence = str(record.seq)
                sequence_id = record.id
                print(f"Processing sequence {sequence_id}...")
                if "N" in original_sequence:
                    continue
                original_nucleotide = original_sequence[position]
                if NO_MUTATION:
                    output_file.write(f">{sequence_id}__{position+1}__{original_nucleotide}to{new_nucleotide}\n")
                    output_file.write(mutated_sequence + "\n")
                else:
                    new_nucleotide = random_mutation(original_nucleotide)
                    mutated_sequence = mutate_sequence(original_sequence, position, new_nucleotide)
                    output_file.write(f">{sequence_id}__{position+1}__{original_nucleotide}to{new_nucleotide}\n")
                    output_file.write(mutated_sequence + "\n")
                counter += 1
                if counter >= COUNTER:
                    break

if __name__ == "__main__":
    os.makedirs(output_dir, exist_ok=True)
    fasta_file = f"{project_root}input_{target}/{type}_{target}.shuffle.fa"  # Replace with your FASTA file path
    main(fasta_file)