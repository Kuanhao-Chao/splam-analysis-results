from Bio import SeqIO
import random
import os
import sys
from pathlib import Path
import pandas as pd

# /ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark
SEQ_LEN = 800
project_root = '/ccb/cybertron/khchao/splam-analysis-results'
test_junc_out_f = f'{project_root}/train/results/tool_benchmark/{SEQ_LEN}bp/d_a.bed'
def main():
    # Positive-MANE
    pos_mane_f = f'{project_root}/train/results/MANE/BAM_REF_Intersection/{SEQ_LEN}bp/100_juncs/d_a.bed'

    # Positive-Alts
    pos_alts_f = f'{project_root}/train/results/RefSeq_ALTS/BAM_REF_Intersection/{SEQ_LEN}bp/100_juncs/d_a.bed'

    # Negative-1
    neg_1_f = f'{project_root}/train/results/Neg_1/Select_junctions/{SEQ_LEN}bp/1_juncs/d_a.bed'

    print("pos_mane_f: ", pos_mane_f)
    print("pos_alts_f: ", pos_alts_f)
    print("neg_1_f: ", neg_1_f)

    files = [pos_mane_f, pos_alts_f, neg_1_f]
    nums = [1000, 1000, 2000]

    # Initialize an empty DataFrame to store the merged data
    merged_df = pd.DataFrame()

    for idx, f in enumerate(files):
        print("f: ", f)
        # Read the BED file
        df = pd.read_csv(f, sep='\t', header=None)

        df = df[df[0].isin(['chr1', 'chr9'])]
        df_sampled = df.sample(n=nums[idx], replace=False)
        if idx == 0 or idx == 1:
            df_sampled['label'] = '+'
        elif idx == 2:
            df_sampled['label'] = '-'
        # print("df_sampled: ", df_sampled)

        # Append the sampled DataFrame with labels to the merged DataFrame
        merged_df = pd.concat([merged_df, df_sampled], ignore_index=True)

    print("Merged DataFrame: ", merged_df)
    # Write the sampled entries to a new BED file
    merged_df.to_csv(test_junc_out_f, sep='\t', index=False, header=False)

    # test_junc_f = f'{project_root}/train/results/tool_benchmark/{SEQ_LEN}bp/d_a.bed'


if __name__ == "__main__":
    main()