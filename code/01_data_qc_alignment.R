# =============================================================================
# PART 1 (Day 1): Data import, QC/EDA, alignment, distance matrix and heatmap
# Output: results/part1_objects.rds  (used by 02_trees_spike.R)
# Run from the repository root (or from code/ - the script handles both).
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
# Install once (uncomment if needed):
# install.packages(c("tidyverse", "rentrez", "ape", "phangorn", "pheatmap"))
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "msa", "DECIPHER"))

library(tidyverse)
library(rentrez)
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(pheatmap)

set.seed(123)

# Always work from the repository root, even if the script is run from code/
if (basename(getwd()) == "code") setwd("..")

dir.create("data/fasta",   showWarnings = FALSE, recursive = TRUE)
dir.create("data/genbank", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",      showWarnings = FALSE, recursive = TRUE)
dir.create("results",      showWarnings = FALSE, recursive = TRUE)
source("code/helpers.R")   # plot_genome_qc()

# ---- 1. Data import ---------------------------------------------------------
meta <- read_csv("data/virus_metadata.csv")
print(meta)

# NCBI allows <= 3 requests per second without an API key
fetch_fasta <- function(acc) {
  f <- file.path("data/fasta", paste0(acc, ".fasta"))
  if (!file.exists(f)) {
    writeLines(entrez_fetch("nuccore", id = acc, rettype = "fasta", retmode = "text"), f)
    Sys.sleep(0.5)
  }
  f
}
fetch_gb <- function(acc) {
  f <- file.path("data/genbank", paste0(acc, ".gb"))
  if (!file.exists(f)) {
    writeLines(entrez_fetch("nuccore", id = acc, rettype = "gb", retmode = "text"), f)
    Sys.sleep(0.5)
  }
  f
}

genomes <- do.call(c, lapply(meta$accession, function(a) readDNAStringSet(fetch_fasta(a))))
names(genomes) <- meta$short_name
print(genomes)

# Download the GenBank records now too (needed on Day 2 for the spike gene)
invisible(lapply(meta$accession, fetch_gb))

# ---- 2. Quality control and exploratory analysis ----------------------------
qc <- meta %>%
  mutate(length  = width(genomes),
         GC      = as.numeric(letterFrequency(genomes, "GC", as.prob = TRUE)) * 100,
         N_bases = as.numeric(letterFrequency(genomes, "N")))
write_csv(qc, "results/genome_qc.csv")
print(qc %>% select(short_name, genus, host, length, GC, N_bases))

stopifnot(all(qc$N_bases < 100))  # sanity check: few ambiguous bases

p_qc <- plot_genome_qc(qc)   # defined in code/helpers.R
print(p_qc)
ggsave("figures/fig1_genome_qc.png", p_qc, width = 8, height = 5, dpi = 300)

# ---- 3. Multiple sequence alignment (nucleocapsid gene) ---------------------
# Aligning whole 30 kb genomes from different genera is very slow and unreliable
# (divergent regions, rearrangements). Instead the nucleocapsid (N) gene, present
# in every coronavirus, is extracted from the GenBank annotation and aligned
# codon-aware (translated, aligned, back-translated) with DECIPHER.
get_cds <- function(acc, genome, pattern) {
  ln <- readLines(fetch_gb(acc))
  starts <- grep("^     CDS ", ln)
  hits <- list()
  for (s in starts) {
    nxt <- which(grepl("^     [A-Za-z]", ln) & seq_along(ln) > s)[1]
    if (is.na(nxt)) nxt <- length(ln) + 1
    block <- paste(trimws(ln[s:(nxt - 1)]), collapse = " ")
    if (grepl(paste0("/product=\"[^\"]*(", pattern, ")"), block, ignore.case = TRUE)) {
      loc  <- sub("^CDS +", "", trimws(ln[s]))
      nums <- as.integer(regmatches(loc, gregexpr("[0-9]+", loc))[[1]])
      if (length(nums) == 2 && !grepl("join|complement", loc))
        hits[[length(hits) + 1]] <- subseq(genome, nums[1], nums[2])
    }
  }
  if (length(hits) == 0) return(NULL)
  hits[[which.max(sapply(hits, length))]]   # keep the longest match
}

n_list <- map2(meta$accession, seq_along(genomes),
               ~ get_cds(.x, genomes[[.y]], "nucleocapsid|nucleoprotein"))
names(n_list) <- meta$short_name
missing <- names(n_list)[map_lgl(n_list, is.null)]
if (length(missing)) message("N gene not found for: ", paste(missing, collapse = ", "))
nset <- DNAStringSet(Filter(Negate(is.null), n_list))
print(width(nset))   # expect roughly 1,100 - 1,500 nt each
meta <- meta[meta$short_name %in% names(nset), ]

aln_file <- "results/N_gene_alignment.fasta"
if (file.exists(aln_file)) {
  aln_set <- readDNAStringSet(aln_file)
} else {
  aln_set <- DECIPHER::AlignTranslation(nset, verbose = FALSE)
  writeXStringSet(aln_set, aln_file)
}

# Trim: remove alignment columns where >50% of sequences have a gap
mat  <- as.matrix(aln_set)
keep <- colMeans(mat == "-") < 0.5
cat("Alignment columns:", ncol(mat), "-> kept after trimming:", sum(keep), "\n")
phy <- phyDat(mat[, keep], type = "DNA")

# ---- 4. Distance matrix and heatmap -----------------------------------------
d_p  <- dist.hamming(phy)          # proportion of differing sites
dmat <- as.matrix(d_p)
stopifnot(!anyNA(dmat))            # heatmap clustering fails if NA/NaN present
write.csv(round(dmat, 4), "results/p_distance_matrix.csv")

annot <- data.frame(Genus = meta$genus, row.names = meta$short_name)
annot <- annot[rownames(dmat), , drop = FALSE]  # same order as the matrix

# Saved directly to file first (most reliable way to save a pheatmap)
pheatmap(dmat, annotation_row = annot,
         main = "Pairwise p-distance of nucleocapsid (N) gene",
         filename = "figures/fig2_distance_heatmap.png", width = 9, height = 7)
# Then draw on screen LAST, so it stays visible in the Plots pane
pheatmap(dmat, annotation_row = annot,
         main = "Pairwise p-distance of nucleocapsid (N) gene")

# ---- 5. Save objects for Part 2 ---------------------------------------------
saveRDS(list(meta = meta, genomes = genomes, phy = phy),
        "results/part1_objects.rds")
cat("Part 1 finished. Objects saved to results/part1_objects.rds\n")
