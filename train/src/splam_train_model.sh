# for seq_len in 200 400 600 800
for seq_len in 40 100
do  
    echo "Processing ${seq_len}bp sequences..."
    python splam_train_model.py $seq_len
done