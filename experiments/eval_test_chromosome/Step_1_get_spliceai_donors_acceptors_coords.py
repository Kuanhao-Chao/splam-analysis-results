from Bio import SeqIO
import random
import os
import sys
from pathlib import Path
import pandas as pd

project_root = '/ccb/cybertron/khchao/splam-analysis-results/'
def main():
    SEQ_LEN = 800
    THRESHOLD = 100
    ref_junc_f = f'{project_root}/train/results/ALL_RefSeq/REF_junctions/ref_d_a.sort.bed'
    test_junc_f = f'{project_root}/train/results/BAM_TEST_junctions/{SEQ_LEN}bp/d_a.sample1.bed'
    junc_fs = [test_junc_f]
    output_dir = f'{project_root}results/eval_test_chromosome/'
    output_files = [output_dir]
    for output_file in output_files:
        os.makedirs(output_file, exist_ok=True)
    nums = [100000]

    COUNTER = 0
    global_df = pd.DataFrame(columns = ['chr', 'start', 'end', 'junc', 'score', 'strand', 'intersect'])
    ref_junc_df = pd.read_csv(ref_junc_f, delimiter="\t", header=None)

    for junc_fidx in range(0, len(junc_fs), 1):
        print("junc_fidx\t: ", junc_fidx)
        junc_f = junc_fs[junc_fidx]
        print("junc_f\t\t: ", junc_f)

        os.makedirs(output_files[junc_fidx]+"splam/", exist_ok=True)
        os.makedirs(output_files[junc_fidx]+"spliceai/", exist_ok=True)

        junc_df = pd.read_csv(junc_f, delimiter="\t", header=None)
        junc_df = junc_df.loc[junc_df[1] > 0]

        print('len(ref_junc_df): ', len(ref_junc_df))
        print('len(junc_df): ', len(junc_df))
        
        # Perform an outer join on the specified columns and add an indicator column
        merged_df = pd.merge(junc_df, ref_junc_df, how='left', on=[0, 1, 2, 5], indicator=True)
        # Add a new column 'is_intersect' to indicate if the entry intersects with ref_junc_df
        merged_df['is_intersect'] = merged_df['_merge'].apply(lambda x: 1 if x == 'both' else 0)
        # Drop the _merge column as it's no longer needed
        merged_df = merged_df.drop(columns=['_merge'])
        # Rename columns to match the desired output format
        out_df = merged_df.rename(columns={0:"chr", 1:"start", 2:"end", "3_x":"junc", "4_x":"score", 5:"strand", 'is_intersect': 'intersect'})
        # Keep only the necessary columns, including the new 'intersect' column
        out_df = out_df[['chr', 'start', 'end', 'junc', 'score', 'strand', 'intersect']]
        print("len(out_df): ", len(out_df))
        print("out_df: ", out_df)
        print(len(out_df[out_df['intersect'] == 1]))

        # if nums[junc_fidx] <= len(junc_df):
        #     junc_df = junc_df.sample(n=nums[junc_fidx], random_state=1).reset_index(drop=True)
        # else:
        #     junc_df = junc_df.sample(n=len(junc_df), random_state=1).reset_index(drop=True)
        #     # junc_df = junc_df.sample(n=nums[junc_fidx], random_state=1, replace=True).reset_index(drop=True)

        ################################
        # SPLAM test data curatiolsn
        ################################
        global_df = out_df
        global_df['end'] -= 1
        print("SPLAM out_df   : ", out_df)
        global_df.to_csv(output_files[junc_fidx]+"splam/splam.juncs.bed", sep="\t", header=None, index=0)

        ################################
        # SpliceAI test data curatiolsn
        ################################
        global_df_spliceai = global_df.copy()
        global_df_spliceai['start'] -= 5200
        global_df_spliceai['end'] += 5200
        print("SpliceAI out_df: ", global_df_spliceai)
        global_df_spliceai.to_csv(output_files[junc_fidx]+"spliceai/spliceai.juncs.bed", sep="\t", header=None, index=0)

if __name__ == "__main__":
    main()