
# =============================================================================
# Helper functions shared by the analysis scripts and the report
# (sourced by 01_data_qc_alignment.R, 02_trees_spike.R and coronavirus_report.Rmd)
# =============================================================================
library(ape)
library(ggplot2)

# ---- Genome length vs GC content plot (Figure 1) ----------------------------
# Uses ggrepel (if installed) so that labels of close points do not overlap.
plot_genome_qc <- function(qc) {
  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(aes(label = short_name), size = 3,
                             show.legend = FALSE, max.overlaps = Inf)
  } else {
    geom_text(aes(label = short_name), vjust = -0.8, size = 3, show.legend = FALSE)
  }
  ggplot(qc, aes(length, GC, colour = genus, shape = host)) +
    geom_point(size = 3) +
    label_layer +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    labs(x = "Genome length (nt)", y = "GC content (%)",
         title = paste("Genome length and GC content of", nrow(qc), "coronaviruses")) +
    theme_minimal()
}

# ---- Bootstrap support for the nodes of a tree ------------------------------
# A split is the division of all tips into two groups made by one branch.
# Written as a text key from the side that does not contain the first tip
# (alphabetically), so the same split gets the same key in every tree.
.canon_split <- function(tips_in_clade, all_tips) {
  side <- if (all_tips[1] %in% tips_in_clade) setdiff(all_tips, tips_in_clade) else tips_in_clade
  paste(sort(side), collapse = "|")
}
.tree_splits <- function(tr, all_tips) {
  pp  <- ape::prop.part(tr)
  lab <- attr(pp, "labels")
  vapply(pp, function(ix) .canon_split(lab[ix], all_tips), character(1))
}

# Returns the bootstrap support (%) of every internal node of `tree`
# (rooted or unrooted), in node order (node Ntip+1 first).
# `bs` is a list of bootstrap trees (class multiPhylo).
# Trivial splits (one tip against the rest, e.g. the outgroup branch) are NA.
node_support <- function(tree, bs) {
  all_tips <- sort(tree$tip.label)
  tree_sp  <- .tree_splits(tree, all_tips)
  bs_sp    <- lapply(seq_along(bs), function(i) .tree_splits(ape::unroot(bs[[i]]), all_tips))
  sup <- vapply(tree_sp, function(s) {
    if (s == "") return(NA_real_)
    100 * mean(vapply(bs_sp, function(b) s %in% b, logical(1)))
  }, numeric(1))
  n_side <- lengths(strsplit(tree_sp, "|", fixed = TRUE))
  sup[n_side < 2] <- NA_real_
  round(sup)
}

# ---- Draw a tree with bootstrap support on the nodes ------------------------
# Only support values >= min_support are printed.
plot_tree_support <- function(tree, bs, tip_col = "black", main = "",
                              min_support = 50, legend_cols = NULL) {
  sup <- node_support(tree, bs)
  op  <- par(mar = c(if (is.null(legend_cols)) 2 else 4, 1, 3, 1))
  on.exit(par(op))
  plot(tree, tip.color = tip_col, main = main, cex = 0.9, label.offset = 0.003)
  lab <- ifelse(!is.na(sup) & sup >= min_support, as.character(sup), "")
  ape::nodelabels(lab, frame = "none", adj = c(1.2, -0.5), cex = 0.7, col = "darkblue")
  ape::add.scale.bar()
  if (!is.null(legend_cols)) {
    legend("bottom", legend = names(legend_cols), text.col = legend_cols,
           horiz = TRUE, bty = "n", xpd = NA, inset = c(0, -0.1), cex = 0.8)
  }
  invisible(sup)
}
