for SPLICEAI_VERSION in 1 2 3 4 5
do
    echo ./2_spliceai_prediction_all_seq_N.sh $SPLICEAI_VERSION &
    ./2_spliceai_prediction_all_seq_N.sh $SPLICEAI_VERSION &
    echo -e "\n"
done
# done