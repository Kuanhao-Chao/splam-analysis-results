# Reviewer Comment
While the architecture of Splam is rooted in convolutional neural networks, recent advancements in deep learning have introduced several tools like EnsembleSplice, Splice2Deep, and SpliceBERT, which employ more sophisticated techniques. These models, such as the bidirectional encoder representations from transformers in SpliceBERT (https://doi.org/10.1101/2023.01.31.526427), offer alternative approaches to solving similar problems. Presenting them alongside Splam provides readers with a comprehensive comparison rather than solely focusing on Splice Ai.

# Tools to explore
- EnsembleSplice (https://github.com/OluwadareLab/EnsembleSplice) -> it seems like the entire pipeline is more geared towards model training than prediction? 
- Splice2Deep (https://github.com/SomayahAlbaradei/Splice_Deep) -> only has pretrained models and works specifically for the given datasets (there are 5 but does not cover mPanTro), not with the adapted inputs
- SpliceBERT (https://github.com/biomed-AI/SpliceBERT?tab=readme-ov-file#how-to-use-splicebert)
- Spliceator (https://git.unistra.fr/nscalzitti/spliceator)
- SpliceFinder (https://github.com/deepomicslab/SpliceFinder)

- Splam
- SpliceAI 

# Datasets
** now linked in my dir
Positive_MANE: /ccb/cybertron/khchao/splam-analysis-results/train/results/MANE/INPUTS/800bp
Positive_ALTS: /ccb/cybertron/khchao/splam-analysis-results/train/results/ALL_RefSeq/INPUTS/800bp
Neg_1: /ccb/cybertron/khchao/splam-analysis-results/train/results/Neg_1/INPUTS/800bp
Neg_Random: /ccb/cybertron/khchao/splam-analysis-results/train/results/Neg_Random/INPUTS/800bp

# Tests

/ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark/spliceai/spliceai.noN.juncs.seq.fa -> location of all dataset
/ccb/cybertron/smao10/splam-analysis-results/train/results/tool_benchmark/spliceai/spliceai.noN.juncs.seq.fa
Write scripts to run Spliceator, SpliceFinder, Splice2Deep, and EnsembleSplice (do in this order), taking a FASTA file as input and outputting a CSV file with junction coordinates, strand, and donor and acceptor scores.


- Standard performance on human dataset
- Generalization test on other species?

Sample 10k from each dataset (alr shuffled) and then run tools
0.5 score threshold -> for positive, negative 

GRAPH
ROC/PR curve 

TABLE
precision/recall/f1/aupr/auroc/top-k accuracy

# setup notes
https://support.idre.ucla.edu/helpdesk/KB/View/67402712-guide-for-installing-tensorflow-and-pytorch -> how i was able to install gpu stuff
https://saturncloud.io/blog/how-to-troubleshoot-pytorchs-torchcudaisavailable-returning-false-in-windows-10/#:~:text=cause%20of%20torch.-,cuda.,able%20to%20detect%20your%20GPU. -> helpful guide
df -H -> disk filesys spaces
du -h -> disk usage