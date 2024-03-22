from Bio import SeqIO

threshold_up = 0.5
threshold_down = 0.01

for threshold in range(2, 101):
    junction_scores = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{threshold}/LOG/_junction_score.bed'

    fasta_file = f'/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_{threshold}/INPUTS/800bp/input_neg_{threshold}_shuffled.fa'

    fasta_outfile_up = f'/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_{threshold}/INPUTS/800bp/input_neg_{threshold}_{threshold_up}_up.fa'

    fasta_outfile_down= f'/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_{threshold}/INPUTS/800bp/input_neg_{threshold}_{threshold_down}_down.fa'

    fasta_outf_up = open(fasta_outfile_up, "w")
    fasta_outf_down = open(fasta_outfile_down, "w")

    # Load FASTA file into a dictionary with SeqIO
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # print("sequences: ", sequences.keys())

    with open(junction_scores, "r") as fr:
        lines = fr.read().split("\n")
        for line in lines:
            eles = line.split("\t")
            if len(eles) < 8:
                continue
            chrname, start, end, junction, _, strand, donor_score, acceptor_score = eles
            
            if float(donor_score) > threshold_up and float(acceptor_score) > threshold_up:
                seq_name = ""
                if strand == "+":
                    seq_name = f"{chrname}:{int(start)-200}-{int(start)+200}({strand})_{chrname}:{int(end)-200}-{int(end)+200}({strand})"
                elif strand == "-":
                    seq_name = f"{chrname}:{int(end)-200}-{int(end)+200}({strand})_{chrname}:{int(start)-200}-{int(start)+200}({strand})"
                # Check if the generated seq_name matches any header in the FASTA dictionary
                if seq_name in sequences.keys():
                    print(f">{seq_name}")
                    print(sequences[seq_name].seq)
                    fasta_outf_up.write(f">{seq_name}\n{sequences[seq_name].seq}\n")

            if float(donor_score) < threshold_down and float(acceptor_score) > threshold_down:
                seq_name = ""
                if strand == "+":
                    seq_name = f"{chrname}:{int(start)-200}-{int(start)+200}({strand})_{chrname}:{int(end)-200}-{int(end)+200}({strand})"
                elif strand == "-":
                    seq_name = f"{chrname}:{int(end)-200}-{int(end)+200}({strand})_{chrname}:{int(start)-200}-{int(start)+200}({strand})"
                # Check if the generated seq_name matches any header in the FASTA dictionary
                if seq_name in sequences.keys():
                    print(f">{seq_name}")
                    print(sequences[seq_name].seq)
                    fasta_outf_down.write(f">{seq_name}\n{sequences[seq_name].seq}\n")

    fasta_outf_up.close()
    fasta_outf_down.close()
