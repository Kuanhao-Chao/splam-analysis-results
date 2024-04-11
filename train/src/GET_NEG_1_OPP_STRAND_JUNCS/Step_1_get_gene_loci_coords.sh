for threshold in {33..100}
do
    JUNC_BEFORE_SORT=/ccb/cybertron/khchao/splam-analysis-results/train/data/junctions_${threshold}.bed
    JUNC=/ccb/cybertron/khchao/splam-analysis-results/train/data/junctions_${threshold}.sort.bed

    sort -k1,1 -k2,2n -k3,3n ${JUNC_BEFORE_SORT} > ${JUNC}

    MANE_dir=/ccb/cybertron/khchao/splam-analysis-results/train/results/MANE/REF_junctions/ref_d_a.sort.bed
    output_dir=/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_${threshold}/BAM_junctions/

    mkdir -p $output_dir

    awk '{if ($6 == "-") {print}}' ${MANE_dir} > ${output_dir}MANE_neg.bed
    awk '{if ($6 == "+") {print}}' ${MANE_dir} > ${output_dir}MANE_pos.bed

    # Finding intervals that are only present in the positive strand file.
    bedtools intersect -sorted -v -a ${output_dir}MANE_pos.bed -b ${output_dir}MANE_neg.bed > ${output_dir}MANE_pos_only.bed

    sort -k1,1 -k2,2n -k3,3n ${output_dir}MANE_pos_only.bed > ${output_dir}MANE_pos_only.sort.bed

    # Finding intervals that are only present in the negative strand file.
    bedtools intersect -sorted -v -a ${output_dir}MANE_neg.bed -b ${output_dir}MANE_pos.bed > ${output_dir}MANE_neg_only.bed

    sort -k1,1 -k2,2n -k3,3n ${output_dir}MANE_neg_only.bed > ${output_dir}MANE_neg_only.sort.bed

    bedtools intersect -sorted -wa -wb -a ${output_dir}MANE_pos_only.bed -b ${JUNC} -sorted -filenames > ${output_dir}neg_hits.bed

    bedtools intersect -sorted -wa -wb -a ${output_dir}MANE_neg_only.bed -b ${JUNC} -sorted -filenames > ${output_dir}pos_hits.bed
done