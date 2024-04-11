SEQ_LEN=$1
project_root=/ccb/cybertron/khchao/splam-analysis-results/
output_dir=${project_root}train/results/ALL_RefSeq/BAM_REF_Intersection/${SEQ_LEN}bp/100_juncs/

bedtools getfasta -s -fi ${project_root}Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ${output_dir}donor.bed -fo ${output_dir}donor_seq.fa

bedtools getfasta -s -fi ${project_root}Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ${output_dir}acceptor.bed -fo ${output_dir}acceptor_seq.fa