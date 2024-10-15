#!/usr/bin/env Rscript
## Script name: save_probe_genes_heatmap.R
## Purpose: source of global scripts and generate output files
## Date Created: October 1, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign input arguments to variables
RUN_ID <- args[1]
PACT_ID <- args[2]
TSV_PATH <- args[3]
METRICS_DIR <- args[4]
PDF_OUT <- args[5]

# Print the input arguments for confirmation
message("RUN_ID:", RUN_ID)
message("PACT_ID:", PACT_ID)
message("TSV_PATH:", TSV_PATH, "\n")
message("METRICS_DIR:", METRICS_DIR, "\n")


# Load necessary libraries ----------------------------------------------------
suppressMessages({
    library("ComplexHeatmap")
    library("tidyverse")
    library("data.table")
    library("grid")
    library("stringr")
    library("circlize")
})


# Function to read coverage files ---------------------------------------------
get_coverage_files <- function() {
    TXT_FILES <- dir(METRICS_DIR, pattern = "PerTargetCoverage.txt",
                     all.files = TRUE, full.names = TRUE)
    message("Reading all txt files...")
    data_list <- lapply(TXT_FILES, function(x) fread(x, data.table = FALSE))
    names(data_list) <- sub('_PerTargetCoverage.txt', '', basename(TXT_FILES))
    return(data_list)
}

# Function to extract gene names ----------------------------------------------
extract_gene_name <- function(row_names) {
    row_names <- sub("_KRAS\\(3845\\)", "_3845(KRAS)", row_names)
    row_names <- sub("PCDHAC2", "(PCDHAC2)", row_names)
    matches <- stringr::str_match(row_names, "\\(([^)]+)\\)")
    genes <- matches[, 2]
    return(genes)
}

# Function to process coverage data -------------------------------------------
get_coverage_data <- function(input_genes) {
    data_list <- get_coverage_files()
    coverage_dfs <- lapply(names(data_list), function(sam_name) {
        df <- data_list[[sam_name]][, c("name", "mean_coverage")]
        colnames(df)[2] <- sam_name
        return(df)
    })

    coverages_df <- base::Reduce(function(x, y)
        base::merge(x, y, by = "name", all = TRUE), coverage_dfs)

    rownames(coverages_df) <- coverages_df$name
    coverages_df$name <- NULL
    new_cols <- stringr::str_split_fixed(colnames(coverages_df), "_", 6)[, 6]
    colnames(coverages_df) <- new_cols
    coverages_df$gene <- extract_gene_name(rownames(coverages_df))
    coverages_df <- coverages_df[coverages_df$gene %in% input_genes &
                                     !is.na(coverages_df$gene), ]
    # Remove Hapmap and NTC columns
    cols2drop <- !grepl("NC_HAPMAP|NTC_H20", colnames(coverages_df))
    coverages_df <- coverages_df[, cols2drop]
    return(coverages_df)
}

# Function grab number of probe rows per gene ---------------------------------
get_gene_chunks <- function(coverage_df) {
    genes_list <- unique(coverage_df$gene)
    total_genes <- length(genes_list)
    genes_per_hm <- ceiling(total_genes / 6)
    return(split(genes_list, ceiling(seq_along(genes_list) / genes_per_hm)))
}

# Function draw individual heatmap page ---------------------------------------
draw_heatmap_page <- function(gene_matrix, hm_title, df_gene) {
    color_scale <- c("navyblue", "white", "darkred")
    color_breaks <- c(0, 200, 500)
    hm_color <- circlize::colorRamp2(color_breaks, color_scale)

    hm <- ComplexHeatmap::Heatmap(
        gene_matrix,
        col = hm_color,
        name = "Mean Coverage",
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_title = hm_title,
        column_title_gp = gpar(fontsize = 16),
        row_split = df_gene$gene,
        row_title_rot = 0,
        heatmap_legend_param = list(
            legend_gp = gpar(fontsize = 12),
            legend_height = unit(6, "cm"),
            at = color_breaks
        )
    )
    return(ComplexHeatmap::draw(hm, padding = unit(c(10, 10, 10, 10), "mm")))
}


# Function to save heatmaps by gene -------------------------------------------
save_heatmap_by_gene <- function(coverage_df, pdf_output, PACT_ID) {
    gene_chunks <- get_gene_chunks(coverage_df)
    hm_title <- paste(PACT_ID, "HS Metrics Select Gene Probe Mean Coverage")

    message("Saving PDF: ", pdf_output)
    pdf(pdf_output, width = 16, height = 8)

    for (genes in gene_chunks) {
        df_gene <- coverage_df[coverage_df$gene %in% genes, ]
        gene_matrix <- as.matrix(df_gene[, !colnames(df_gene) %in% "gene"])
        gene_matrix <- gene_matrix[complete.cases(gene_matrix), ]

        if (nrow(gene_matrix) == 0) next

        draw_heatmap_page(gene_matrix, hm_title, df_gene)
    }

    dev.off()
}

# Function to generate heatmaps for runs --------------------------------------
generate_run_heatmaps <- function(PACT_ID, RUN_ID, TSV_PATH, PDF_OUT) {
    input_genes <- read.delim(TSV_PATH, sep = "\t",
                              header = F, stringsAsFactors = F)[[1]]
    message(paste("Running:", PACT_ID, RUN_ID))
    coverage_df <- get_coverage_data(input_genes)
    pdf_output <- PDF_OUT
    #pdf_output <- paste(PACT_ID, RUN_ID, "Mean_Coverage_probe_hm.pdf", sep = "_")
    save_heatmap_by_gene(coverage_df, pdf_output, PACT_ID)
    message("Heatmap Generation Complete!")
}

# Main execution of the script ------------------------------------------------
generate_run_heatmaps(PACT_ID, RUN_ID, TSV_PATH, PDF_OUT)
