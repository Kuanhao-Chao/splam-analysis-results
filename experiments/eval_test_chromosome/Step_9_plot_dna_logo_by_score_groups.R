library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# Neg_Xs
spliceai_idx <- "3"
for (tool in c('splam', 'spliceai')) {
# for (tool in c('spliceai')) {
  for (target in c("up", "down")) {
    for (ref_label in c("0", "1")) {

      if (tool == 'splam') {
        fasta_file <- paste("/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/", tool, "/predict_out/", tool ,"/score_", target, "_", ref_label, ".fa", sep = "")
        output_dir <- paste("/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/", tool, "/predict_out/", tool, sep = "")
      } else if (tool == 'spliceai') {
        fasta_file <- paste("/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/", tool, "/predict_out/", tool , spliceai_idx, "/score_", target, "_", ref_label, ".fa", sep = "")
        output_dir <- paste("/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/", tool, "/predict_out/", tool, spliceai_idx, sep = "")
      }
      print(fasta_file)
      dna <- readDNAStringSet(fasta_file)
      sample_num <- min(50000, length(dna))  # Use the smaller of 50000 or the number of sequences
      print(output_dir)
      # Check if output directory exists before creating
      if (!dir.exists(output_dir)) {
        dir.create(output_dir)
      }

      if (tool == 'splam') {
        if (sample_num > 0 && all(width(dna[1:sample_num]) >= 620)) {  # Check that we have sequences and they are long enough
          # Donor
          ggseqlogo(as.character(subseq(dna[1:sample_num], 580, 620)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
          ggsave(
            paste(output_dir, "/acceptor", target, "_", ref_label, ".png", sep = ""),  # Note: This is labeled as "acceptor" but is generated from the "donor" region
            width = 15,
            height = 3,
            dpi = 300
          )
          # Acceptor
          ggseqlogo(as.character(subseq(dna[1:sample_num], 180, 220)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
          ggsave(
            paste(output_dir, "/donor", target, "_", ref_label, ".png", sep = ""),  # Note: This is labeled as "donor" but is generated from the "acceptor" region
            width = 15,
            height = 3,
            dpi = 300
          )
        } else {
          warning(paste("Skipping", fasta_file, "due to insufficient number of sequences or sequence length."))
        }
      } else if (tool == 'spliceai') {
          # Assuming `dna` is a list or vector of DNA sequences
          # Preprocess to extract the last 5220 to 5180 region for each sequence
          preprocessed_dna_acceptor <- lapply(dna, function(seq) {
            seq_length <- nchar(seq) # Get the length of the sequence
            if (seq_length >= 5220) {
              start_idx <- seq_length - 5220 + 1
              end_idx <- seq_length - 5180  
              # print(length(substring(seq, start_idx, end_idx)))
              return(as.character(substring(seq, start_idx, end_idx)))
            } else {
              return(NULL) # Exclude sequences that are too short
            }
          })
          # Filter out NULL values if any sequences were too short
          preprocessed_dna_acceptor <- Filter(Negate(is.null), preprocessed_dna_acceptor)
          # Flatten the list to a character vector if necessary
          preprocessed_dna_acceptor <- unlist(preprocessed_dna_acceptor)
          # Filter sequences to keep only those with exactly 400 characters
          print((preprocessed_dna_acceptor))
          # Check if we have any sequences left after filtering
          if (length(preprocessed_dna_acceptor) > 0) {
            # Plot the DNA logo for the acceptor site
            dna_logo_acceptor <- ggseqlogo(preprocessed_dna_acceptor, method = 'bits', seq_type = 'dna', col_scheme = 'nucleotide')
            print(dna_logo_acceptor)
            # Save the plot
            ggsave(paste(output_dir, "/acceptor", target, "_", ref_label, ".png", sep = ""),
                  plot = dna_logo_acceptor,
                  width = 15,
                  height = 3,
                  dpi = 300)
          } else {
            warning("No sequences were long enough for the specified region.")
          }

          # Donor
          start_idx = 5180
          end_idx = 5220
  # print(as.character(subseq(dna[1:sample_num], start_idx, end_idx)))
          ggseqlogo(as.character(subseq(dna[1:sample_num], start_idx, end_idx)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
          ggsave(
            paste(output_dir, "/donor", target, "_", ref_label, ".png", sep = ""),  # Note: This is labeled as "donor" but is generated from the "acceptor" region
            width = 15,
            height = 3,
            dpi = 300
          )
        } else {
          warning(paste("Skipping", fasta_file, "due to insufficient number of sequences or sequence length."))
        }
    }
  }
}
