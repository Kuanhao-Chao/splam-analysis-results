for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287 ; do
    echo scp /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/ribozero/$sample/junction_score.bed s3:/ccb/salz3/kh.chao/SPLAM/results/ribozero/$sample/
    scp /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/ribozero/$sample/junction_score.bed s3:/ccb/salz3/kh.chao/SPLAM/results/ribozero/$sample/
done
