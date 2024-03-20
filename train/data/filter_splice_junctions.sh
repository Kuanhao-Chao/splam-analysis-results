for threshold in {23..99}
do
    awk -v bound="$threshold" '{if($5 == bound && $6 != "?") {print}}' ./all.def.junctions.bed > ./junctions_${threshold}.bed
done