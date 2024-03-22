library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# Neg_Xs
for (x in 2:100) {
  for (target in c("up", "down")) {
    if (target == "up") {
      fasta_file <- paste("/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_", x, "/INPUTS/800bp/input_neg_", x, "_0.5_up.fa", sep = "")
    } else {
      fasta_file <- paste("/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_", x, "/INPUTS/800bp/input_neg_", x, "_0.01_down.fa", sep = "")
    }

    print(fasta_file)
    dna <- readDNAStringSet(fasta_file)

    sample_num <- min(50000, length(dna))  # Use the smaller of 50000 or the number of sequences
    
    if (target == "up") {
      output_dir <- paste("negs/neg_", x, "/", target, "_0.5", sep = "")
    } else if (target == "down") {
      output_dir <- paste("negs/neg_", x, "/", target, "_0.01", sep = "")
    }
    print(output_dir)

    # Check if output directory exists before creating
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }

    if (sample_num > 0 && all(width(dna[1:sample_num]) >= 620)) {  # Check that we have sequences and they are long enough
      # Donor
      ggseqlogo(as.character(subseq(dna[1:sample_num], 580, 620)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
      ggsave(
        paste(output_dir, "/acceptor.png", sep = ""),  # Note: This is labeled as "acceptor" but is generated from the "donor" region
        width = 15,
        height = 3,
        dpi = 300
      )
      # Acceptor
      ggseqlogo(as.character(subseq(dna[1:sample_num], 180, 220)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
      ggsave(
        paste(output_dir, "/donor.png", sep = ""),  # Note: This is labeled as "donor" but is generated from the "acceptor" region
        width = 15,
        height = 3,
        dpi = 300
      )
    } else {
      warning(paste("Skipping", fasta_file, "due to insufficient number of sequences or sequence length."))
    }
  }
}
