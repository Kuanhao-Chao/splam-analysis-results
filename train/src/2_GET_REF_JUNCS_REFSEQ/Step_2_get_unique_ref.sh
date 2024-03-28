extracted_junc_type="all"
if [ "$extracted_junc_type" == "all" ]; then
    project_root=/ccb/cybertron/khchao/splam-analysis-results/
    output_dir=${project_root}train/results/ALL_RefSeq/REF_junctions/
    sort -k1,1 -k2,2n -k3,3n ${output_dir}ref_d_a_${extracted_junc_type}.bed | uniq > ${output_dir}ref_d_a_${extracted_junc_type}.sort.bed
else
    project_root=/ccb/cybertron/khchao/splam-analysis-results/
    output_dir=${project_root}train/results/ALL_RefSeq/REF_junctions/
    sort -k1,1 -k2,2n -k3,3n ${output_dir}ref_d_a.bed | uniq > ${output_dir}ref_d_a.sort.bed
fi
