import os
def main():
    for threshold in range(2, 101):
        output_dir=f'/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_{threshold}/BAM_junctions/'
        input_files = [f'{output_dir}neg_hits.bed', f'{output_dir}pos_hits.bed']
        ofile = f'{output_dir}junctions_{threshold}_cleaned.bed'
        with open(ofile, "w") as fw:
            for input_file_name in input_files:
                input_file_path = os.path.join(output_dir, input_file_name)
                print("file: ", input_file_path)
                with open(input_file_path, "r") as f:
                    for line in f:
                        columns = line.strip().split("\t")
                        # print(columns)
                        gtf_start, gtf_end = int(columns[1]), int(columns[2])
                        junc_start, junc_end = int(columns[7]), int(columns[8])
                        gtf_strand, junc_strand = columns[5], columns[11]

                        if gtf_start < junc_start and gtf_end > junc_end and gtf_strand != junc_strand:
                            fw.write("\t".join(columns[6:]) + "\n")

if __name__ == "__main__":
    main()