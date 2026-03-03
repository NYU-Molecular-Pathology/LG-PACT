#!/usr/bin/env Rscript
## Script name: exome_cvg_200_csv.R
## Purpose: Intakes run_id to read PerTargetCoverage and saves >200 rc
## Date Created: February 20, 2026
## Version: 1.1.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

args <- commandArgs(TRUE)
RUN_ID      <- args[1]  # e.g. "260211_NB551709_0365_AH2LTGBGY2"
PACT_ID     <- args[2]  # e.g. "PACT-26-6"
METRICS_DIR <- args[3]  # path to CollectHsMetrics directory
DEMUX_CSV   <- args[4]  # path to demux-samplesheet.csv
SAMPLES_CSV <- args[5]  # path to samples.tumor.normal.csv

suppressMessages({
    library("dplyr")
    library("readr")
    library("stringr")
})

# Output to current working directory (Nextflow handles publishDir)
OUTPUT_DIR <- "."

# Functions -------------------------------------------------------------------

# Returns all the per-target coverage files from the staged metrics directory
get_hs_files <- function(metrics_dir) {
    list.files(path = metrics_dir,
               pattern = "PerTargetCoverage\\.txt$",
               full.names = TRUE)
}

# Reads PerTargetCoverage files & saves matrix using file name for SampleID
get_probe_matrix <- function(hs_files) {
    all_rows <- lapply(hs_files, function(f) {
        raw_name  <- str_remove(basename(f), "_PerTargetCoverage\\.txt$")
        sampleID  <- str_split_fixed(raw_name, "_", 6)[1, 6]
        sampleID  <- stringr::str_replace_all(sampleID, "AAH3HVMM5_", "")
        readr::read_tsv(f, col_types = cols(), show_col_types = FALSE) %>%
            dplyr::filter(read_count < 200) %>%
            dplyr::mutate(SampleID = sampleID) %>%
            dplyr::select(SampleID, chrom, start, end, length, name, read_count)
    })
    all_rows <- bind_rows(all_rows)

    probeMat <- all_rows %>% arrange(SampleID) %>%
        filter(SampleID != "NTC", !str_detect(chrom, "chrY"))
    return(probeMat)
}

# Returns the pairs of tumor and normal samples from the staged CSV
get_sample_pairs <- function(samples_csv) {
    all_pairs <- read.csv(samples_csv)
    all_pairs$Normal <- str_split_fixed(all_pairs$Normal, "_", 6)[, 6]
    all_pairs$Tumor  <- str_split_fixed(all_pairs$Tumor,  "_", 6)[, 6]
    all_pairs$Tumor  <- stringr::str_replace_all(all_pairs$Tumor,  "AAH3HVMM5_", "")
    all_pairs$Normal <- stringr::str_replace_all(all_pairs$Normal, "AAH3HVMM5_", "")
    return(all_pairs)
}

# Separate the matrix into two data frames for Normal and Tumor samples
extract_pair_dfs <- function(mat, pairs_df) {
    rows_normal <- str_detect(mat$SampleID, paste0(pairs_df$Normal, collapse = "|"))
    rows_tumor  <- str_detect(mat$SampleID, paste0(pairs_df$Tumor,  collapse = "|"))
    list(normal = as.data.frame(mat[rows_normal, ]),
         tumor  = as.data.frame(mat[rows_tumor,  ]))
}

# Save the matrix to a CSV file
# NOTE: file names must match the Nextflow output block exactly:
#   {RUN_ID}_Normal_PerTargetCoverage_200.csv
#   {RUN_ID}_Tumor_PerTargetCoverage_200.csv
save_matrix <- function(run_id, perTargCovg_probes, sample_type) {
    output_path <- file.path(
        OUTPUT_DIR,
        paste(run_id, sample_type, "PerTargetCoverage_200.csv", sep = "_")
    )
    message("Saving file: ", output_path)
    readr::write_csv(perTargCovg_probes, output_path)
}

# Divides the total probes by the number of probes with reads below 200
# NOTE: file names must match the Nextflow output block exactly:
#   {pact_id}_<normals/tumors>_percentages_200.csv
calc_above_200 <- function(df, pact_id, total_probes, sam_type) {
    freq_table <- table(df$SampleID)
    counts_df <- as.data.frame(freq_table, stringsAsFactors = FALSE)
    counts_df$PACT_ID <- pact_id
    counts_df$Freq <- total_probes - as.numeric(counts_df$Freq)
    counts_df$Freq <- counts_df$Freq / total_probes * 100
    counts_df$Freq <- round(counts_df$Freq, 2)
    colnames(counts_df) <- c("SampleID", "PercentCoverage", "Run")

    csv_out <- paste(pact_id, sam_type, "percentages_200.csv", sep = "_")
    csv_path <- file.path(OUTPUT_DIR, csv_out)
    message(paste0(capture.output(counts_df), collapse = "\n"))
    message("Saving file: ", csv_path)
    write.csv(counts_df, file = csv_path, row.names = FALSE, quote = FALSE)
}

# Main execution --------------------------------------------------------------
run_hs_files       <- get_hs_files(METRICS_DIR)
sample_pairs_truth <- get_sample_pairs(SAMPLES_CSV)
perTargCovg_truth  <- get_probe_matrix(run_hs_files)

both_runs <- list(
    truth = extract_pair_dfs(perTargCovg_truth, sample_pairs_truth)
)

save_matrix(RUN_ID, both_runs$truth$normal, "Normal")
save_matrix(RUN_ID, both_runs$truth$tumor,  "Tumor")

# Get total probes minus last six rows (chrY)
total_probes <- nrow(readr::read_tsv(run_hs_files[1],
                                     show_col_types = FALSE)) - 6

calc_above_200(both_runs$truth$normal, PACT_ID, total_probes, "normals")
calc_above_200(both_runs$truth$tumor,  PACT_ID, total_probes, "tumors")