#!/usr/bin/env Rscript

library(argparse)
library(ape)
library(glue)
library(dplyr)
library(grid)
library(ComplexHeatmap)
#library(signals)
library(magick)

# Suppress ComplexHeatmap messages
ht_opt$message = FALSE
library(here)
source(here("src/heatmap_plot.R"))

#' Generate a Heatmap Aligned with a Phylogenetic Tree
#' 
#' @param treefile Path to the tree file (Newick format)
#' @param clusters Optional data frame with cell_id column for subsetting
#' @param cnv_data Data frame containing CNV data with cell_id column (required only if clusters not provided)
#' @param output_file Optional path to save the output plot
#' @param plot_width Width of the plot (default: 89 * 0.039)
#' @param plot_height Height of the plot (default: 2)
#' @param chroms Character vector of chromosome labels
#' @param title Plot title (default: empty string)
#' @param tree_width Width of the tree portion (default: 1)
#' @param linkheight Height of the linking lines (default: 2)
#' @param clone_colors Named vector of colors for clones
#' @return List containing heatmap object and tree object
make_heatmap_tree <- function(treefile, 
                              clusters = NULL, 
                              cnv_data = NULL,
                              output_file = NULL,
                              plot_width = 89 * 0.039,
                              plot_height = 2,
                              chroms = c(paste0(1:11), "13", "15", "17", "20", "X"),
                              title = "",
                              tree_width = 1,
                              linkheight = 2,
                              clone_colors = c("A" = "firebrick4", "B" = "deepskyblue4")) {
  
  # Validate inputs
  if (!file.exists(treefile)) {
    stop("Tree file does not exist: ", treefile)
  }
  
  # Check if cnv_data is required
  if (is.null(clusters) && is.null(cnv_data)) {
    stop("Either clusters or cnv_data must be provided")
  }
  
  # Read the tree
  cat("Reading tree from:", treefile, "\n")
  mytree <- ape::read.tree(file = treefile)
  cat("Tree has", length(mytree$tip.label), "tips\n")
  
  # Remove "cell_" prefix from tree tip labels if present
  if (any(grepl("^cell_", mytree$tip.label))) {
    cat("Removing 'cell_' prefix from tree tip labels\n")
    mytree$tip.label <- gsub("^cell_", "", mytree$tip.label)
    cat("Updated tree tip labels (first 3):", paste(head(mytree$tip.label, 3), collapse = ", "), "\n")
  }
  
  # Handle clusters and CNV data logic
  if (!is.null(clusters)) {
    if (!"cell_id" %in% colnames(clusters)) {
      stop("clusters data frame must contain 'cell_id' column")
    }
    cat("Filtering tree to", nrow(clusters), "cells from clusters\n")
    mytree <- ape::keep.tip(mytree, clusters$cell_id)
    
    # If clusters is provided but cnv_data is not, try to use global dat$cn
    if (is.null(cnv_data)) {
      if (exists("dat") && "cn" %in% names(dat)) {
        cnv_data <- dat$cn
      } else {
        stop("When clusters is provided without cnv_data, global dat$cn must exist")
      }
    }
  }
  
  # Validate cnv_data
  if (!"cell_id" %in% colnames(cnv_data)) {
    stop("cnv_data must contain 'cell_id' column")
  }
  
  cat("CNV data has", nrow(cnv_data), "rows and", length(unique(cnv_data$cell_id)), "unique cells\n")
  
  # Check overlap between tree and CNV data
  overlap_cells <- intersect(mytree$tip.label, unique(cnv_data$cell_id))
  cat("Found", length(overlap_cells), "overlapping cells between tree and CNV data\n")
  
  if (length(overlap_cells) == 0) {
    cat("Tree tips (first 5):", paste(head(mytree$tip.label, 5), collapse = ", "), "\n")
    cat("CNV cell_ids (first 5):", paste(head(unique(cnv_data$cell_id), 5), collapse = ", "), "\n")
    stop("No overlapping cells found between tree and CNV data")
  }
  
  # Filter the CNV data to match the cells in the tree
  cnaneuploid <- cnv_data %>% 
    filter(cell_id %in% mytree$tip.label)
  
  cat("Filtered CNV data has", nrow(cnaneuploid), "rows for", length(unique(cnaneuploid$cell_id)), "cells\n")
  
  # Fix mixed data types that cause plotting errors
  cnaneuploid <- cnaneuploid %>%
    mutate(
      # Ensure consistent character types
      chr = as.character(chr),
      cell_id = as.character(cell_id),
      patient = as.character(patient),
      Data_Type = as.character(Data_Type),
      Sample_Name = as.character(Sample_Name),
      
      # Handle columns that might have mixed types
      Genotype = as.character(Genotype),
      Group = as.character(Group),
      
      # Ensure numeric columns are properly numeric
      state = as.numeric(as.character(state)),
      start = as.numeric(as.character(start)),
      end = as.numeric(as.character(end)),
      segment_size_bp = as.numeric(as.character(segment_size_bp)),
      segment_size_MB = as.numeric(as.character(segment_size_MB))
    )
  
  cat("Data types standardized successfully\n")
  
  # NEW: Chromosome diagnostic
  cat("=== CHROMOSOME DIAGNOSTIC ===\n")
  actual_chroms <- sort(unique(cnaneuploid$chr))
  expected_chroms <- chroms
  
  cat("Chromosomes in data:", paste(actual_chroms, collapse = ", "), "\n")
  cat("Expected chromosomes:", paste(expected_chroms, collapse = ", "), "\n")
  
  missing_in_data <- setdiff(expected_chroms, actual_chroms)
  extra_in_data <- setdiff(actual_chroms, expected_chroms)
  
  if (length(missing_in_data) > 0) {
    cat("Missing from data:", paste(missing_in_data, collapse = ", "), "\n")
  }
  if (length(extra_in_data) > 0) {
    cat("Extra in data:", paste(extra_in_data, collapse = ", "), "\n")
  }
  
  # Auto-adjust chroms to match data
  if (length(missing_in_data) > 0 || length(extra_in_data) > 0) {
    cat("🔧 Auto-adjusting chromosome list to match data\n")
    chroms <- actual_chroms
    cat("Updated chroms:", paste(chroms, collapse = ", "), "\n")
  }
  cat("============================\n")
  
  cat("Chromosome column type:", class(cnaneuploid$chr), "\n")
  
  if (nrow(cnaneuploid) == 0) {
    stop("No CNV data remains after filtering to match tree cells")
  }
  
  # Check if plotHeatmap function exists
  if (!exists("plotHeatmap")) {
    stop("plotHeatmap function not found. Please ensure it's loaded in your environment.")
  }
  
  # --- Plotting ---
  p <- plotHeatmap(
    cnaneuploid,
    tree = mytree %>% ape::ladderize(),
    column_title = title,
    column_title_gp = gpar(fontsize = 8),
    linkheight = linkheight,
    chrlabels = chroms,
    show_heatmap_legend = FALSE,
    plotfrequency = FALSE,
    frequency_height = 0.5,
    anno_width = 0.02,
    annofontsize = 7,
    show_legend = FALSE,
    show_clone_text = FALSE,
    show_library_label = FALSE,
    show_clone_label = FALSE,
    plottree = TRUE,
    reorderclusters = TRUE,
    tree_width = tree_width,
    clone_pal = clone_colors,
    clusters = clusters
  )
  
  # Capture the plot as a grid object
  pout <- grid.grabExpr(draw(p), width = plot_width, height = plot_height)
  
  # Save output if requested
  if (!is.null(output_file)) {
    # Determine file type from extension
    ext <- tools::file_ext(output_file)
    
    if (ext %in% c("pdf", "PDF")) {
      pdf(output_file, width = plot_width, height = plot_height)
      draw(p)
      dev.off()
    } else if (ext %in% c("png", "PNG")) {
      png(output_file, width = plot_width * 100, height = plot_height * 100, res = 300)
      draw(p)
      dev.off()
    } else {
      warning("Unsupported file format. Supported formats: pdf, png")
    }
  }
  
  return(list(hm = pout, tree = mytree))
}

# Main function with argparse
main <- function() {
  # Set up for both interactive and command-line use
  if (interactive()) {
    # Mimic command-line input for testing
    argv <- c("-t", "tree.nwk",
              "-c", "clusters.csv",
              "-o", "output.pdf",
              "--width", "5",
              "--height", "3",
              "--title", "Test Heatmap")
  } else {
    # Get real command-line arguments
    argv <- commandArgs(trailingOnly = TRUE)
  }
  
  # Create argument parser
  parser <- ArgumentParser(description = "Generate a heatmap aligned with a phylogenetic tree")
  
  # Required arguments
  parser$add_argument("-t", "--treefile", type = "character", required = TRUE,
                      help = "Path to the tree file (Newick format)")
  
  # Optional arguments
  parser$add_argument("-c", "--clusters", type = "character", default = NULL,
                      help = "CSV file with cell_id column for subsetting")
  
  parser$add_argument("-d", "--cnv_data", type = "character", default = NULL,
                      help = "CSV file with CNV data (required only if clusters not provided)")
  
  parser$add_argument("-o", "--output", type = "character", default = NULL,
                      help = "Output file path (pdf or png)")
  
  parser$add_argument("--width", type = "double", default = 89 * 0.039,
                      help = "Plot width (default: 3.471)")
  
  parser$add_argument("--height", type = "double", default = 2,
                      help = "Plot height (default: 2)")
  
  parser$add_argument("--title", type = "character", default = "",
                      help = "Plot title (default: empty)")
  
  parser$add_argument("--tree_width", type = "double", default = 1,
                      help = "Width of the tree portion (default: 1)")
  
  parser$add_argument("--linkheight", type = "double", default = 2,
                      help = "Height of the linking lines (default: 2)")
  
  parser$add_argument("--chroms", type = "character", nargs = "*",
                      default = NULL,
                      help = "Chromosome labels (space-separated). If not provided, will auto-detect from data")
  
  # Parse arguments
  args <- parser$parse_args(argv)
  
  # Load data files
  clusters <- NULL
  cnv_data <- NULL
  
  if (!is.null(args$clusters)) {
    if (!file.exists(args$clusters)) {
      stop("Clusters file does not exist: ", args$clusters)
    }
    clusters <- read.csv(args$clusters)
  }
  
  if (!is.null(args$cnv_data)) {
    if (!file.exists(args$cnv_data)) {
      stop("CNV data file does not exist: ", args$cnv_data)
    }
    cnv_data <- read.csv(args$cnv_data)
  }
  
  # Validate that either clusters or cnv_data is provided
  if (is.null(clusters) && is.null(cnv_data)) {
    stop("Either --clusters or --cnv_data must be provided")
  }
  
  # Auto-detect or use provided chromosomes
  if (is.null(args$chroms) && !is.null(cnv_data)) {
    # Auto-detect chromosomes from data
    chroms_to_use <- sort(unique(as.character(cnv_data$chr)))
    cat("Auto-detected chromosomes from data:", paste(chroms_to_use, collapse = ", "), "\n")
  } else if (!is.null(args$chroms)) {
    chroms_to_use <- args$chroms
  } else {
    # Fallback to default
    chroms_to_use <- c(paste0(1:11), "13", "15", "17", "20", "X")
  }
  
  # Set up clone colors
  clone_colors <- c("A" = "firebrick4", "B" = "deepskyblue4")
  
  # Call the function
  tryCatch({
    result <- make_heatmap_tree(
      treefile = args$treefile,
      clusters = clusters,
      cnv_data = cnv_data,
      output_file = args$output,
      plot_width = args$width,
      plot_height = args$height,
      chroms = chroms_to_use,  # Fixed: use chroms_to_use instead of args$chroms
      title = args$title,
      tree_width = args$tree_width,
      linkheight = args$linkheight,
      clone_colors = clone_colors
    )
    
    cat("Heatmap tree generated successfully!\n")
    if (!is.null(args$output)) {
      cat("Output saved to:", args$output, "\n")
    }
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}
# Run main function if script is executed directly
if (!interactive()) {
  main()
}