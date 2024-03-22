library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# # Mane only
# dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_MANE.shuffle.fa")
# dir.create("MANE")
# # Donor
# ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
# ggsave(
#   "MANE/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
# ggsave(
#   "MANE/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )


# # Refseq alternative only
# dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_ALTS.shuffle.fa")
# dir.create("ALTS")
# # Donor
# ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
# ggsave(
#   "ALTS/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
# ggsave(
#   "ALTS/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )


library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# Neg_Xs
for (x in 2:100) {
  # Refseq alternative only
  fasta_file <- paste("/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_", x, "/INPUTS/800bp/input_neg_", x, "_shuffled.fa", sep = "")
  print(fasta_file)
  dna <- readDNAStringSet(fasta_file)

  sample_num <- min(50000, length(dna))  # Use the smaller of 50000 or the number of sequences
  
  output_dir <- paste("negs/neg_", x, sep = "")
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



# # Refseq alternative only
# dna <- readDNAStringSet("../../../src/INPUTS/800bp/input_neg_random.shuffle.fa")
# dir.create("neg_random")
# # Donor
# print("as.character(subseq(dna[1:50000], 580, 620)): ", as.character(subseq(dna[1:50000], 580, 620)))
# ggseqlogo(as.character(subseq(dna[1:50000], 580, 620)), method = 'probability', seq_type='dna', col_scheme='nucleotide')
# ggseqlogo(as.character(subseq(dna[1:50000], 580, 620)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
# ggsave(
#   "neg_random/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna[1:50000], 180, 220)), method = 'probability', seq_type='dna', col_scheme='nucleotide')
# ggseqlogo(as.character(subseq(dna[1:50000], 180, 220)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
# ggsave(
#   "neg_random/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
