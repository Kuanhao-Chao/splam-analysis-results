project_root=/ccb/cybertron/khchao/splam-analysis-results/
output_dir=${project_root}train/results/ALL_RefSeq/REF_junctions/
sort -k1,1 -k2,2n -k3,3n ${output_dir}ref_d_a.bed | uniq > ${output_dir}ref_d_a.sort.bed