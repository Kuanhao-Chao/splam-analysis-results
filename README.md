# splam-analysis-results
all the scripts to reproduce the results in the splam paper


## Generalization Tests:
`benchmark/src_generalization_test/`

- all outputs from every step of the data processing will be saved in the corresponding folder of the step `{#}_output/`
- `{name}` refers to the name of the gene database

## Generate and pre-process the data:

1. Positive Dataset

    The positive data is extracted from the complete genomic GFF annotation files of each species. 

    1. To start, retrieve all introns from the GFF file:
            
            $ python 1_get_pos_set.py

        Inputs:
        - `{name}.gff` = annotations corresponding to a genome
        - `{name}_genomic.fa` = genome's fasta file 
        - `{name}_annotation_report.txt` = annotation report 
        *These were downloaded together from the NCBI Genome Database*

        Outputs:
        - `databases/{name}.db` = sqlite3-style databases parsed from the `.gff` annotation files
        - `{name}_introns.bed` = protein-coding genes
    
    2. Extract all the sequences from the introns, with the tailored specifications of Splam. This step also performs checks and filters to ensure high quality positive data.

            $ python 2_extract_pos_splam.py
        
        Inputs:
        - `{name}_introns.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `donor.bed`, `acceptor.bed` = bed files referring to the specific 400nt donor and acceptor sequences
        - `donor_seq.fa`, `acceptor_seq.fa`= fasta files containing the 400nt sequences
        - `d_a.bed` = bed file referring to the intron coordinates of the splice junction (midpoint of the donor and acceptor coords)
        - `input_neg_random.fa` = fasta file containing the 800nt sequence that is given to Splam

    3. Extract the same sequences, but with the specifications for SpliceAI.

            $ python 3_extract_pos_spliceai.py

        Inputs:
        - `d_a.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `coords.bed` = bed file referring to the start and end positions of the *whole* SpliceAI input
        - `seq_noN.fa` = fasta file containing the SpliceAI input, with the full flanking sequence (coords refer to splice junction)
        - `seq_N.fa` = fasta file containing the SpliceAI input, with repeating N flanking sequence (coords refer to splice junction)




2. Negative Dataset
    
    You will need to randomly generate the negative splice junction dataset. This works by taking the existing protein-coding genes, then selecting the opposite strand to guarantee unique sequences, and creating pseudo-splice-junctions from GT-AG pairs found on this strand.
    
    1. To start, retrieve the protein-coding genes from all four genomes, making use of the `gffutils` library: 

            $ python 1_get_neg_set.py

        Inputs:
        - `{name}.gff` = annotations corresponding to a genome
        - `{name}_genomic.fa` = genome's fasta file 
        - `{name}_annotation_report.txt` = annotation report 
        *These were downloaded together from the NCBI Genome Database*

        Outputs:
        - `databases/{name}.db` = sqlite3-style databases parsed from the `.gff` annotation files
        - `{name}_genes.bed` = protein-coding genes

    2. Then, generate the dataset of splice junctions, and process into a format readable by Splam. 
        
            $ python 2_extract_neg_splam.py

        Inputs:
        - `{name}_genes.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `donor.bed`, `acceptor.bed` = bed files referring to the specific 400nt donor and acceptor sequences
        - `donor_seq.fa`, `acceptor_seq.fa`= fasta files containing the 400nt sequences
        - `d_a.bed` = bed file referring to the intron coordinates of the splice junction (midpoint of the donor and acceptor coords)
        - `input_neg_random.fa` = fasta file containing the 800nt sequence that is given to Splam

    3. Now process the same splice junctions into a format readable by SpliceAI (it will take a subset of the junctions for efficiency).

            $ python 3_extract_neg_spliceai.py
        
        Inputs:
        - `d_a.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `coords.bed` = bed file referring to the start and end positions of the *whole* SpliceAI input
        - `seq_noN.fa` = fasta file containing the SpliceAI input, with the full flanking sequence (coords refer to splice junction)
        - `seq_N.fa` = fasta file containing the SpliceAI input, with repeating N flanking sequence (coords refer to splice junction)


## Run Splam and SpliceAI on the processed transcripts, and record the scores

Run the following three steps in both folders of the pipeline. They are essentially the same for both positive and negative datasets.

1. Run Splam.

        $ ./4_splam_runner.sh
    
    Inputs:
    - `input_neg_random.fa` 

    Outputs:
    - `score.bed` = bed file containing the Splam-scored splice junctions


2. Run SpliceAI. Depending on your system, this may take several days, so you can run each dataset separately: 

        $ ./5_spliceai_prediction_wrapper.sh {name}
    
    Inputs:
    - `seq_noN.fa`
    - `seq_N.fa`

    Outputs:
    *There are 5 model output folders, each containing 4 folders with the database names*
    - `spliceai_all_seq.name.noN.{name}.tsv`, `spliceai_all_seq.name.N.{name}.tsv` = names and identifiers for the scored splice junctions
    - `spliceai_all_seq.score.a.noN.{name}.tsv`, `spliceai_all_seq.score.a.N.{name}.tsv` = acceptor site scores for every nt in sequence
    - `spliceai_all_seq.score.d.noN.{name}.tsv`, `spliceai_all_seq.score.d.N.{name}.tsv` = donor site scores for every nt in sequence
    - `spliceai_all_seq.score.n.noN.{name}.tsv`, `spliceai_all_seq.score.n.N.{name}.tsv` = neutral (neither) scores for every nt in sequence

3. Post-process the Splam and SpliceAI scores into a single file for comparison:

        $ 6_compile_data.py
    
    Inputs:
    - `seq_noN.fa`, `seq_N.fa`
    - `spliceai_all_seq.name.noN.{name}.tsv`, `spliceai_all_seq.name.N.{name}.tsv`
    - `spliceai_all_seq.score.a.noN.{name}.tsv`, `spliceai_all_seq.score.a.N.{name}.tsv`
    - `spliceai_all_seq.score.d.noN.{name}.tsv`, `spliceai_all_seq.score.d.N.{name}.tsv`
    - `score.bed`

    Outputs:
    - `Splam/{name}.splam_data.csv` = the compiled Splam scores with various identifiers and metrics
    - `SpliceAI/spliceai_data.v{#}.noN.{name}.csv`, `SpliceAI/spliceai_data.v{#}.noN.{name}.csv` = the compiled SpliceAI scores with various identifiers and metrics
    - `combine/aggregate_data.noN.{name}.csv`, `combine/aggregate_data.N.{name}.csv` = the complete list of data compiled from both Splam and SpliceAI, with all data included (no averaging)
    - `combine/averaged_data.noN.{name}.csv`, `combine/averaged_data.N.{name}.csv` = averaged list across the 5 models, most interfaceable form


## Plotting Results:

Navigate to the `figures/` folder. Here, you have a subdirectory for every figure we generated. Simply navigate to the figure you want, and run the `plot.py` function.

* For the `f1_scores`, run the `stats.py` which will generate a `results.csv` for various performance metrics at a specified threshold (that you can edit in the code).

