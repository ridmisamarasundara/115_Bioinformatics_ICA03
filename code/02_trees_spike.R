# =============================================================================
# PART 2 (Day 2): NJ/UPGMA trees, ML tree, spike gene tree, session info
# Requires: results/part1_objects.rds (created by 01_data_qc_alignment.R)
# Run from the repository root (or from code/ - the script handles both).
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
library(tidyverse)
library(Biostrings)
library(msa)
library(ape)
library(phangorn)

set.seed(123)
n_boot <- 50   # number of bootstrap replicates (raise to 100 if time allows)

if (basename(getwd()) == "code") setwd("..")
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
source("code/helpers.R")   # node_support() and plot_tree_support()

# Helper: show a base-R plot on screen AND save it as a PNG in figures/
plot_and_save <- function(file, plot_fun, w = 8, h = 6) {
  plot_fun()
  png(file.path("figures", file), width = w, height = h, units = "in", res = 300)
  plot_fun()
  dev.off()
}

# Load Part 1 results
obj     <- readRDS("results/part1_objects.rds")
meta    <- obj$meta
genomes <- obj$genomes
phy     <- obj$phy

# ---- 1. Distance-based trees (NJ, UPGMA) with bootstrap ---------------------
dm         <- dist.ml(phy, model = "JC69")
nj_tree    <- NJ(dm)
upgma_tree <- upgma(dm)

outgroup <- "IBV"  # Gammacoronavirus = outgroup for alpha + beta
nj_root  <- root(nj_tree, outgroup = outgroup, resolve.root = TRUE)

bs_nj <- bootstrap.phyDat(phy, function(x) NJ(dist.ml(x, model = "JC69")), bs = n_boot)

col_genus <- c(Alphacoronavirus = "darkgreen", Betacoronavirus = "firebrick",
               Gammacoronavirus = "steelblue")
tip_col <- col_genus[meta$genus[match(nj_root$tip.label, meta$short_name)]]

plot_and_save("fig3_nj_tree.png", function() {
  plot_tree_support(nj_root, bs_nj, tip_col = tip_col, legend_cols = col_genus,
                    main = paste0("Neighbour-joining tree, N gene (", n_boot, " bootstraps, JC69)"))
})

# UPGMA assumes a constant rate - shown for comparison only
plot(upgma_tree, main = "UPGMA tree (assumes a constant rate)")
png("figures/fig6_upgma_tree.png", width = 8, height = 6, units = "in", res = 300)
plot(upgma_tree, main = "UPGMA tree (assumes a constant rate)")
dev.off()

# ---- 2. Maximum-likelihood tree ---------------------------------------------
# Compare substitution models using BIC (restricted set for speed)
mt <- modelTest(phy, tree = nj_tree, model = c("JC", "HKY", "GTR"), multicore = FALSE)
mt <- mt[order(mt$BIC), ]
write_csv(mt, "results/model_test.csv")
print(head(mt))
best_model <- mt$Model[1]
print(best_model)

base_model <- sub("\\+.*", "", best_model)  # e.g. "GTR" from "GTR+G+I"
fit0 <- update(pml(nj_tree, phy), k = 4, inv = 0.2)
fit  <- optim.pml(fit0, model = base_model, optGamma = TRUE, optInv = TRUE,
                  rearrangement = "NNI", control = pml.control(trace = 0))
print(fit$logLik)

# Bootstrap (number of replicates set by n_boot at the top)
bs_ml    <- bootstrap.pml(fit, bs = n_boot, optNni = TRUE, multicore = FALSE,
                          control = pml.control(trace = 0))
ml_root  <- root(fit$tree, outgroup = outgroup, resolve.root = TRUE)
tip_col2 <- col_genus[meta$genus[match(ml_root$tip.label, meta$short_name)]]

plot_and_save("fig4_ml_tree.png", function() {
  plot_tree_support(ml_root, bs_ml, tip_col = tip_col2, legend_cols = col_genus,
                    main = paste0("Maximum-likelihood tree, N gene (", best_model, ", ", n_boot, " bootstraps)"))
})
write.tree(ml_root, "results/ml_tree.nwk")

# Do NJ and ML agree? Robinson-Foulds distance (0 = identical topology)
print(RF.dist(unroot(nj_root), unroot(ml_root)))

# ---- 3. Focused analysis: spike (S) gene ------------------------------------
# Extract the spike CDS from the cached GenBank annotation (data/genbank/)
get_spike <- function(acc, genome) {
  ln <- readLines(file.path("data/genbank", paste0(acc, ".gb")))
  starts <- grep("^     CDS ", ln)
  for (s in starts) {
    nxt <- which(grepl("^     [A-Za-z]", ln) & seq_along(ln) > s)[1]
    if (is.na(nxt)) nxt <- length(ln) + 1
    block <- paste(trimws(ln[s:(nxt - 1)]), collapse = " ")
    if (grepl("/product=\"[^\"]*(spike|surface glycoprotein)", block, ignore.case = TRUE)) {
      loc  <- sub("^CDS +", "", trimws(ln[s]))
      nums <- as.integer(regmatches(loc, gregexpr("[0-9]+", loc))[[1]])
      if (length(nums) == 2 && !grepl("join|complement", loc))
        return(subseq(genome, nums[1], nums[2]))
    }
  }
  NULL
}

spike_list <- map2(meta$accession, seq_along(genomes),
                   ~ get_spike(.x, genomes[[.y]]))
names(spike_list) <- meta$short_name
missing <- names(spike_list)[map_lgl(spike_list, is.null)]
if (length(missing)) message("Spike not found for: ", paste(missing, collapse = ", "))
spike <- DNAStringSet(Filter(Negate(is.null), spike_list))
print(width(spike))

# Spike is highly divergent between genera: use beta-coronaviruses only
beta    <- meta$short_name[meta$genus == "Betacoronavirus"]
spike_b <- spike[names(spike) %in% beta]

aln_s <- msa(spike_b, method = "ClustalOmega", type = "DNA")
mat_s <- as.matrix(as(aln_s, "DNAStringSet"))
phy_s <- phyDat(mat_s[, colMeans(mat_s == "-") < 0.5], type = "DNA")

nj_s <- NJ(dist.ml(phy_s, model = "JC69"))
bs_s <- bootstrap.phyDat(phy_s, function(x) NJ(dist.ml(x, model = "JC69")), bs = n_boot)
nj_s <- root(nj_s, outgroup = "HCoV-OC43", resolve.root = TRUE)  # embecovirus outgroup

plot_and_save("fig5_spike_tree.png", function() {
  plot_tree_support(nj_s, bs_s,
                    main = paste0("Spike gene NJ tree (Betacoronavirus, ", n_boot, " bootstraps)"))
})

# ---- 4. Save objects for the report (coronavirus_report.Rmd) ----------------
saveRDS(list(nj_root = nj_root, bs_nj = bs_nj, upgma_tree = upgma_tree,
             mt = mt, best_model = best_model, ml_root = ml_root, bs_ml = bs_ml,
             logLik = fit$logLik, rf = RF.dist(unroot(nj_root), unroot(ml_root)),
             nj_s = nj_s, bs_s = bs_s, spike_width = width(spike),
             spike_missing = missing, n_boot = n_boot),
        "results/part2_objects.rds")

# ---- 5. Session information -------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
sessionInfo()
