#!/usr/bin/env Rscript
# Dengue Variant Tracker - Quality Control and Motif Analysis
# Processes dengue virus sequences and identifies mutation patterns

# Load required libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
  library(dplyr)
  library(ggplot2)
})

cat("==========================================\n")
cat("Dengue Variant Tracker - QC & Analysis\n")
cat("==========================================\n\n")

# Configuration
USE_TEST_DATA <- TRUE  # Set to FALSE for full dataset (requires more RAM)
INPUT_FILE <- if (USE_TEST_DATA) "data/raw/test_dengue.fasta" else "data/raw/dengue_sequences.fasta"
OUTPUT_DIR <- "data/processed"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
if (!dir.exists("plots")) {
  dir.create("plots", recursive = TRUE)
}

# ===== STEP 1: Load Sequences =====
cat("[Step 1/5] Loading sequences...\n")
if (!file.exists(INPUT_FILE)) {
  stop(paste("ERROR: Input file not found:", INPUT_FILE, 
             "\nPlease run ./download_data.sh first"))
}

seqs <- readDNAStringSet(INPUT_FILE)
cat(sprintf("✓ Loaded %d sequences\n", length(seqs)))
cat(sprintf("  Total nucleotides: %s\n", format(sum(width(seqs)), big.mark=",")))

# ===== STEP 2: Quality Control =====
cat("\n[Step 2/5] Performing quality control...\n")

# Calculate basic statistics
seq_widths <- width(seqs)
gc_content <- letterFrequency(seqs, "GC", as.prob = TRUE)
n_content <- letterFrequency(seqs, "N", as.prob = TRUE)

# Quality filtering criteria
MIN_LENGTH <- 500
MAX_N_PERCENT <- 0.05  # 5% maximum ambiguous bases

# Apply filters
cat("  Filtering criteria:\n")
cat(sprintf("    - Minimum length: %d bp\n", MIN_LENGTH))
cat(sprintf("    - Maximum N content: %.1f%%\n", MAX_N_PERCENT * 100))

filtered_seqs <- seqs[seq_widths >= MIN_LENGTH & n_content <= MAX_N_PERCENT]

cat(sprintf("✓ Quality filtering complete\n"))
cat(sprintf("  Sequences retained: %d/%d (%.1f%%)\n", 
            length(filtered_seqs), length(seqs),
            100 * length(filtered_seqs) / length(seqs)))

# Generate QC summary
qc_summary <- data.frame(
  Sequence_ID = names(filtered_seqs),
  Length = width(filtered_seqs),
  GC_Content = as.numeric(letterFrequency(filtered_seqs, "GC", as.prob = TRUE)),
  N_Content = as.numeric(letterFrequency(filtered_seqs, "N", as.prob = TRUE)),
  A_Count = as.numeric(letterFrequency(filtered_seqs, "A")),
  T_Count = as.numeric(letterFrequency(filtered_seqs, "T")),
  G_Count = as.numeric(letterFrequency(filtered_seqs, "G")),
  C_Count = as.numeric(letterFrequency(filtered_seqs, "C"))
)

write.csv(qc_summary, file.path(OUTPUT_DIR, "qc_summary.csv"), row.names = FALSE)
cat(sprintf("✓ QC summary saved: %s\n", file.path(OUTPUT_DIR, "qc_summary.csv")))

# ===== STEP 3: Motif Analysis =====
cat("\n[Step 3/5] Analyzing motifs and mutation patterns...\n")

# Define dengue-specific motifs of interest
# Based on literature on dengue virus mutations and vaccine escape
motifs <- list(
  "ATG" = "Start codon (translation initiation)",
  "GAC" = "Common mutation site in envelope (E) protein",
  "AATAAA" = "Polyadenylation signal",
  "CACAG" = "5' UTR conserved region",
  "AGAGA" = "3' UTR cyclization sequence element"
)

cat("  Scanning for known dengue motifs:\n")
motif_results <- list()

for (motif_seq in names(motifs)) {
  cat(sprintf("    - %s: %s\n", motif_seq, motifs[[motif_seq]]))
  
  # Search for motif in all sequences
  matches <- vmatchPattern(motif_seq, filtered_seqs)
  
  # Convert to data frame
  match_data <- lapply(seq_along(matches), function(i) {
    if (length(matches[[i]]) > 0) {
      data.frame(
        Sequence_ID = names(filtered_seqs)[i],
        Motif = motif_seq,
        Description = motifs[[motif_seq]],
        Position = start(matches[[i]]),
        Count = length(matches[[i]]),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  
  motif_results[[motif_seq]] <- do.call(rbind, match_data)
}

# Combine all motif results
all_motifs <- do.call(rbind, motif_results)
if (is.null(all_motifs)) {
  all_motifs <- data.frame(
    Sequence_ID = character(),
    Motif = character(),
    Description = character(),
    Position = integer(),
    Count = integer()
  )
}

# Calculate motif summary statistics
motif_summary <- all_motifs %>%
  group_by(Motif, Description) %>%
  summarise(
    Total_Occurrences = n(),
    Sequences_With_Motif = n_distinct(Sequence_ID),
    Avg_Occurrences_Per_Seq = mean(Count),
    .groups = 'drop'
  ) %>%
  arrange(desc(Total_Occurrences))

write.csv(all_motifs, file.path(OUTPUT_DIR, "motif_matches.csv"), row.names = FALSE)
write.csv(motif_summary, file.path(OUTPUT_DIR, "motif_summary.csv"), row.names = FALSE)

cat(sprintf("✓ Motif analysis complete\n"))
cat(sprintf("  Total motif matches found: %d\n", nrow(all_motifs)))
cat(sprintf("  Results saved: %s\n", file.path(OUTPUT_DIR, "motif_matches.csv")))

# ===== STEP 4: Generate Visualizations =====
cat("\n[Step 4/5] Generating visualizations...\n")

# Plot 1: Sequence length distribution
p1 <- ggplot(qc_summary, aes(x = Length)) +
  geom_histogram(bins = 30, fill = "#2E86AB", color = "white", alpha = 0.8) +
  geom_vline(xintercept = MIN_LENGTH, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Dengue Virus Sequence Length Distribution",
    subtitle = sprintf("n = %d sequences (after QC filtering)", nrow(qc_summary)),
    x = "Sequence Length (bp)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave("plots/01_sequence_length_dist.png", p1, width = 8, height = 5, dpi = 300)
cat("  ✓ Saved: plots/01_sequence_length_dist.png\n")

# Plot 2: GC content distribution
p2 <- ggplot(qc_summary, aes(x = GC_Content * 100)) +
  geom_histogram(bins = 20, fill = "#A23B72", color = "white", alpha = 0.8) +
  geom_density(aes(y = after_stat(count)), color = "#F18F01", size = 1.2) +
  labs(
    title = "GC Content Distribution",
    subtitle = "Dengue virus genome composition",
    x = "GC Content (%)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave("plots/02_gc_content_dist.png", p2, width = 8, height = 5, dpi = 300)
cat("  ✓ Saved: plots/02_gc_content_dist.png\n")

# Plot 3: Motif occurrence frequency
if (nrow(motif_summary) > 0) {
  p3 <- ggplot(motif_summary, aes(x = reorder(Motif, Total_Occurrences), 
                                   y = Total_Occurrences)) +
    geom_col(fill = "#06A77D", alpha = 0.8) +
    geom_text(aes(label = Total_Occurrences), hjust = -0.2, size = 4) +
    coord_flip() +
    labs(
      title = "Motif Frequency in Dengue Sequences",
      subtitle = "Known mutation and regulatory sites",
      x = "Motif",
      y = "Total Occurrences"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(family = "mono")
    )
  
  ggsave("plots/03_motif_frequency.png", p3, width = 8, height = 5, dpi = 300)
  cat("  ✓ Saved: plots/03_motif_frequency.png\n")
}

# Plot 4: Nucleotide composition
nucleotide_comp <- qc_summary %>%
  summarise(
    A = sum(A_Count),
    T = sum(T_Count),
    G = sum(G_Count),
    C = sum(C_Count)
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "Nucleotide", values_to = "Count")

p4 <- ggplot(nucleotide_comp, aes(x = Nucleotide, y = Count, fill = Nucleotide)) +
  geom_col(alpha = 0.8, color = "white", size = 1) +
  scale_fill_manual(values = c(A = "#E63946", T = "#457B9D", G = "#F4A261", C = "#2A9D8F")) +
  labs(
    title = "Overall Nucleotide Composition",
    subtitle = "Across all dengue sequences",
    x = "Nucleotide",
    y = "Total Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave("plots/04_nucleotide_composition.png", p4, width = 7, height = 5, dpi = 300)
cat("  ✓ Saved: plots/04_nucleotide_composition.png\n")

# ===== STEP 5: Save processed sequences =====
cat("\n[Step 5/5] Saving processed data...\n")

writeXStringSet(filtered_seqs, file.path(OUTPUT_DIR, "filtered_sequences.fasta"))
cat(sprintf("✓ Filtered sequences saved: %s\n", 
            file.path(OUTPUT_DIR, "filtered_sequences.fasta")))

# Create analysis summary
analysis_summary <- list(
  timestamp = Sys.time(),
  input_file = INPUT_FILE,
  total_sequences_loaded = length(seqs),
  sequences_after_qc = length(filtered_seqs),
  qc_filter_rate = round(100 * length(filtered_seqs) / length(seqs), 1),
  min_length_threshold = MIN_LENGTH,
  max_n_threshold = MAX_N_PERCENT,
  mean_sequence_length = round(mean(seq_widths), 0),
  mean_gc_content = round(mean(gc_content) * 100, 2),
  total_motif_matches = nrow(all_motifs),
  unique_motifs_found = length(unique(all_motifs$Motif))
)

# Save as JSON-like text for easy reading
sink(file.path(OUTPUT_DIR, "analysis_summary.txt"))
cat("DENGUE VARIANT TRACKER - ANALYSIS SUMMARY\n")
cat("==========================================\n\n")
for (name in names(analysis_summary)) {
  cat(sprintf("%-30s: %s\n", name, analysis_summary[[name]]))
}
sink()

# Memory cleanup
gc()

cat("\n==========================================\n")
cat("Analysis Complete!\n")
cat("==========================================\n")
cat("\nOutput files created:\n")
cat("  - data/processed/qc_summary.csv\n")
cat("  - data/processed/motif_matches.csv\n")
cat("  - data/processed/motif_summary.csv\n")
cat("  - data/processed/filtered_sequences.fasta\n")
cat("  - data/processed/analysis_summary.txt\n")
cat("  - plots/*.png (4 visualization files)\n")
cat("\nNext step:\n")
cat("  Launch dashboard: Rscript -e \"shiny::runApp('app.R')\"\n")
cat("==========================================\n")
