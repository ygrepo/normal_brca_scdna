#!/usr/bin/env Rscript

library(ape)
library(glue)
library(dplyr)
library(grid)
library(ComplexHeatmap)

#' Generate a Heatmap Aligned with a Phylogenetic Tree
#' 
#' @param treefile Path to the tree file (Newick format)
#' @param clusters Optional data frame with cell_id column for subsetting
#' @param cnv_data Data frame containing CNV data with cell_id column
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
  
  if (is.null(cnv_data)) {
    stop("CNV data must be provided")
  }
  
  # Read the tree
  mytree <- ape::read.tree(file = treefile)
  
  # Optionally, subset the tree to only include cells from provided clusters
  if (!is.null(clusters)) {
    if (!"cell_id" %in% colnames(clusters)) {
      stop("clusters data frame must contain 'cell_id' column")
    }
    mytree <- ape::keep.tip(mytree, clusters$cell_id)
  }
  
  # Get cell order from the tree
  cellorder <- mytree$tip.label[mytree$edge[, 2]]
  
  # Filter the CNV data to match the cells in the tree
  if (!"cell_id" %in% colnames(cnv_data)) {
    stop("cnv_data must contain 'cell_id' column")
  }
  
  cnaneuploid <- cnv_data %>% 
    filter(cell_id %in% mytree$tip.label)
  
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

#' Command line interface for make_heatmap_tree
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript heatmap_tree.R <treefile> [options]\n")
    cat("Options:\n")
    cat("  --clusters <file>     CSV file with cell_id column\n")
    cat("  --cnv_data <file>     CSV file with CNV data (required)\n")
    cat("  --output <file>       Output file (pdf or png)\n")
    cat("  --width <num>         Plot width (default: 3.471)\n")
    cat("  --height <num>        Plot height (default: 2)\n")
    cat("  --title <string>      Plot title\n")
    cat("  --help                Show this help message\n")
    quit(status = 1)
  }
  
  # Initialize parameters
  treefile <- args[1]
  clusters <- NULL
  cnv_data <- NULL
  output_file <- NULL
  plot_width <- 89 * 0.039
  plot_height <- 2
  title <- ""
  
  # Parse optional arguments
  i <- 2
  while (i <= length(args)) {
    if (args[i] == "--clusters") {
      clusters <- read.csv(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--cnv_data") {
      cnv_data <- read.csv(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--output") {
      output_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--width") {
      plot_width <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--height") {
      plot_height <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--title") {
      title <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--help") {
      cat("Usage: Rscript heatmap_tree.R <treefile> [options]\n")
      cat("Options:\n")
      cat("  --clusters <file>     CSV file with cell_id column\n")
      cat("  --cnv_data <file>     CSV file with CNV data (required)\n")
      cat("  --output <file>       Output file (pdf or png)\n")
      cat("  --width <num>         Plot width (default: 3.471)\n")
      cat("  --height <num>        Plot height (default: 2)\n")
      cat("  --title <string>      Plot title\n")
      cat("  --help                Show this help message\n")
      quit(status = 0)
    } else {
      i <- i + 1
    }
  }
  
  # Call the function
  result <- make_heatmap_tree(
    treefile = treefile,
    clusters = clusters,
    cnv_data = cnv_data,
    output_file = output_file,
    plot_width = plot_width,
    plot_height = plot_height,
    title = title
  )
  
  cat("Heatmap tree generated successfully!\n")
  if (!is.null(output_file)) {
    cat("Output saved to:", output_file, "\n")
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}