
echo Generating scores 
mkdir -p ../results/splam/2/

echo python 2_splam_predict.py -f ../results/splam/1/neg_input.fa -o ../results/splam/2/neg_score.bed -m ../../../model/splam_script.pt
python 2_splam_predict.py -f ../results/splam/1/neg_input.fa -o ../results/splam/2/neg_score.bed -m ../../../model/splam_script.pt

echo python 2_splam_predict.py -f ../results/splam/1/pos_input.fa -o ../results/splam/2/pos_score.bed -m ../../../model/splam_script.pt
python 2_splam_predict.py -f ../results/splam/1/pos_input.fa -o ../results/splam/2/pos_score.bed -m ../../../model/splam_script.pt

