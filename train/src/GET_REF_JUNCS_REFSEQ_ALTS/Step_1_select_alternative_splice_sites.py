import os 
import re
import pandas as pd
import sys

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
output_dir = f'{project_root}train/results/RefSeq_ALTS/REF_junctions/'

def chr_name_convert():
    f_chrs = open(f'{project_root}/Dataset/Refseq_2_UCSU_chromosome_names.tsv', "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = eles[1]
    return chrs


def main():
    JUNC_COUNTER = 0
    THRESHOLD = "100"
    SEQ_LEN=sys.argv[1]
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4
    os.makedirs(output_dir, exist_ok=True)
    d_a_out = f'{output_dir}ref_d_a.bed'

    MANE_output_dir = f'{project_root}train/results/MANE/BAM_REF_Intersection/{SEQ_LEN}bp/{THRESHOLD}_juncs/'
    ALL_RefSeq_output_dir = f'{project_root}train/results/ALL_RefSeq/BAM_REF_Intersection/{SEQ_LEN}bp/{THRESHOLD}_juncs/'
    MANE_juncs = pd.read_csv(f'{MANE_output_dir}d_a.bed', sep="\t", header=None)
    refseq_juncs = pd.read_csv(f'{ALL_RefSeq_output_dir}d_a.bed', sep="\t", header=None)
    print(MANE_juncs)
    print(refseq_juncs)
    intersect_df = pd.merge(refseq_juncs, MANE_juncs, how ='left', on =[0, 1, 2, 5])
    intersect_df = intersect_df[intersect_df.isnull().any(axis=1)]
    print(intersect_df)
    # intersect_df=intersect_df.dropna()
    intersect_df = intersect_df.drop(['3_y', '4_y'], axis=1)
    out_df = intersect_df.rename(columns={0:"chr",1:"start", 2:"end", "3_x":"junc", "4_x":"score", 5:"strand"})
    out_df.to_csv(d_a_out, sep="\t", index=False, header=None)
    

if __name__ == "__main__":
    main()