import pandas as pd

################################
# SpliceAI output
################################
from Bio import SeqIO

threshold_up = 0.5
threshold_down = 0.1

scores_csv="/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/merged.csv"
scores_df = pd.read_csv(scores_csv, sep='\t', header=0)
print("scores_df: ", scores_df)

spliceai_model = 3
for target in ['splam', 'spliceai']:
    print("target: ", target)
    for ref_label_target in [0, 1]:
        score_up_count = 0
        score_down_count = 0
        # Load FASTA file into a dictionary with SeqIO
        fasta_file=f'/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/{target}/{target}.juncs.seq.fa'
        print("fasta_file: ", fasta_file)
        if target == 'splam':
            fasta_outfile_up = f'/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/{target}/predict_out/{target}/score_up_{ref_label_target}.fa'
            fasta_outfile_down=f'/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/{target}/predict_out/{target}/score_down_{ref_label_target}.fa'
        elif target == 'spliceai':
            fasta_outfile_up = f'/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/{target}/predict_out/{target}{spliceai_model}/score_up_{ref_label_target}.fa'
            fasta_outfile_down = f'/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/{target}/predict_out/{target}{spliceai_model}/score_down_{ref_label_target}.fa'     
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        fasta_outf_up = open(fasta_outfile_up, "w")
        fasta_outf_down = open(fasta_outfile_down, "w")
        
        with open(scores_csv, "r") as fr:
            lines = fr.read().split("\n")
            for line in lines[1:]:
                eles = line.split("\t")
                if len(eles) < 9:
                    continue
                chrname = eles[0]
                start = eles[1]
                end = eles[2]
                strand = eles[3]
                spliceai_score_a = eles[4]
                spliceai_score_d = eles[5]
                splam_score_a = eles[6]
                splam_score_d = eles[7]
                ref_label = eles[8]
                if target == 'splam':
                    donor_score = float(splam_score_d)
                    acceptor_score = float(splam_score_a)
                elif target == 'spliceai':
                    donor_score = float(spliceai_score_d)
                    acceptor_score = float(spliceai_score_a)
                # print("donor_score: ", donor_score)
                # print("acceptor_score: ", acceptor_score)

                if float(donor_score) > threshold_up and float(acceptor_score) > threshold_up and float(ref_label) == ref_label_target:
                    score_up_count += 1
                    seq_name = ""
                    if target == 'splam':
                        if strand == "+":
                            seq_name = f"{chrname}:{int(start)-200}-{int(start)+200}({strand})_{chrname}:{int(end)-200}-{int(end)+200}({strand})"
                        elif strand == "-":
                            seq_name = f"{chrname}:{int(end)-200}-{int(end)+200}({strand})_{chrname}:{int(start)-200}-{int(start)+200}({strand})"
                    elif target == 'spliceai':
                        seq_name = f"{chrname}:{int(start)-5200}-{int(end)+5200}({strand})"
                        print("seq_name: ", seq_name)
                    # Check if the generated seq_name matches any header in the FASTA dictionary
                    if seq_name in sequences.keys():
                        # print(f">{seq_name}")
                        # print(sequences[seq_name].seq)
                        fasta_outf_up.write(f">{seq_name}\n{sequences[seq_name].seq}\n")

                if float(donor_score) < threshold_down and float(acceptor_score) < threshold_down and float(ref_label) == ref_label_target:
                    score_down_count += 1
                    seq_name = ""
                    if target == 'splam':
                        if strand == "+":
                            seq_name = f"{chrname}:{int(start)-200}-{int(start)+200}({strand})_{chrname}:{int(end)-200}-{int(end)+200}({strand})"
                        elif strand == "-":
                            seq_name = f"{chrname}:{int(end)-200}-{int(end)+200}({strand})_{chrname}:{int(start)-200}-{int(start)+200}({strand})"
                    elif target == 'spliceai':
                        seq_name = f"{chrname}:{int(start)-5200}-{int(end)+5200}({strand})"
                        print("seq_name: ", seq_name)
                    # Check if the generated seq_name matches any header in the FASTA dictionary
                    if seq_name in sequences.keys():
                        # print(f">{seq_name}")
                        # print(sequences[seq_name].seq)
                        fasta_outf_down.write(f">{seq_name}\n{sequences[seq_name].seq}\n")
        fasta_outf_up.close()
        fasta_outf_down.close()

print("score_up_count: ", score_up_count)
print("score_down_count: ", score_down_count)
