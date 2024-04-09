# for OUTPUTFILE in "./dataset/pos" "./dataset/neg_1" "./dataset/neg_random"
for OUTPUTFILE in "/ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark/"
do
    for TARGET in "splam" "spliceai"
    # for TARGET in "spliceai"
    do
        if [[ $TARGET == "splam" ]]
        then
            echo "bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.seq.fa"

            bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.seq.fa

            echo "bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.seq.fa"

            bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.seq.fa
        fi

        if [[ $TARGET == "spliceai" ]]
        then
            echo "bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.seq.fa"

            bedtools getfasta -s -fi /ccb/cybertron/khchao/splam-analysis-results/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.seq.fa
        fi
    done
done