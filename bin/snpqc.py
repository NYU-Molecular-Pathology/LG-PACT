#!/usr/bin/env python3
"""
Generate SNPQC : Creates marker ratio file, correlation matrix heatmap and SNP VAF plot

"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import re

# Function to process individual pileup files
def process_pileup_file(file_path):
    # Read the pileup file into a DataFrame (space-separated, no header)
    try:
        data = pd.read_csv(file_path, sep=" ", header=None)
    except pd.errors.EmptyDataError:
        print(f"The file {file_path} is empty, an empty DataFrame was created.")
        return pd.DataFrame()
    else:
        # Create POS by concatenating the first two columns (chromosome and position)
        data["POS"] = data[0].astype(str) + ":" + data[1].astype(str)
    
        # Column 2 is used as the reference allele (REF)
        data["REF"] = data[2]
    
        # Calculate the REF ratio if read depth > 50; else assign NaN
        data["Ratio"] = data.apply(lambda row: row[3].count(row["REF"]) / len(row[3]) 
                                if isinstance(row[3], str) and len(row[3]) > 50 else np.nan, axis=1)
    
        # Return DataFrame with POS, REF, and calculated ratio
    return data[["POS", "REF", "Ratio"]]

# Function to merge pileup data from multiple files
def merge_dataframes(dir_path,marker_bedfile):
    # Load marker positions from the BED file
    marker_df = pd.read_csv(marker_bedfile, sep="\t", header=None)
    marker_df["POS"] = marker_df[0].astype(str) + ":" + marker_df[2].astype(str)
    markers = marker_df[["POS"]]
    
    final_df = pd.DataFrame()
    
    # Walk through directory and process .pileup files
    for root, dirs, files in os.walk(dir_path):
        for file in files:
            if file.endswith(".pileup"):
                file_path = os.path.join(root, file)
                data = process_pileup_file(file_path)
                if data.empty:
                    continue
                else:
                    data.rename(columns={"Ratio": os.path.splitext(file)[0]}, inplace=True)
                    # Merge processed data
                    if final_df.empty:
                        final_df = data
                    else:
                        final_df = pd.merge(final_df, data, on=["POS", "REF"], how='outer')
    
    # Merge with marker positions to ensure all positions are retained
    final_df = pd.merge(markers, final_df, on=["POS"], how="outer")
    return final_df

# Function to create marker ratio file and save to CSV
def create_marker_ratiofile(path_to_pileup,rundir,marker_bedfile):
    final_dataframe = merge_dataframes(path_to_pileup,marker_bedfile)
    
    # Filter for chromosomes starting with 'chr'
    final_dataframe = final_dataframe[final_dataframe['POS'].astype(str).str.startswith('chr')]
    
    # Separate control and other columns
    control_columns = [col for col in final_dataframe.columns if col.startswith('0_')]
    other_columns = [col for col in final_dataframe.columns if not col.startswith('0_')]
    
    other_columns_s = sorted(other_columns[2:])
    
    # Rearrange columns: first two, remove controls, then the rest
    new_order = other_columns[:2] + other_columns_s
    final_dataframe_reordered = final_dataframe[new_order]
    
    # Save to CSV
    clinical_snpqc_path = os.path.join(rundir,"clinical/snpqc") 
    if not os.path.exists(clinical_snpqc_path): os.mkdir(clinical_snpqc_path,0o0755)

    final_dataframe_reordered.to_csv(os.path.join(clinical_snpqc_path, "phase3_1000_GRCh37_markers_ratio.csv"), index=False, na_rep='NA')
    print("Final dataframe saved under", clinical_snpqc_path)
    return final_dataframe_reordered

# Create heatmap of correlation between samples
def create_correlation_matrix_heatmap(marker_ratio_df, runid, rundir):
    correlation_matrix = np.zeros((marker_ratio_df.shape[1], marker_ratio_df.shape[1]))
    
    for i in range(marker_ratio_df.shape[1]):
        for j in range(marker_ratio_df.shape[1]):
            if i == j:
                correlation_matrix[i][j] = 1
            else:
                calculate_diff = np.abs(marker_ratio_df.iloc[:, i] - marker_ratio_df.iloc[:, j]) < 0.1
                either_nan = np.isnan(marker_ratio_df.iloc[:, i]) | np.isnan(marker_ratio_df.iloc[:, j])
                
                if np.sum(either_nan) > 0.2 * marker_ratio_df.shape[0]:
                    correlation_matrix[i][j] = np.nan
                else:
                    valid_comparison = calculate_diff & ~either_nan
                    same_measure_count = np.sum(valid_comparison)
                    total_positions = marker_ratio_df.shape[0] - np.sum(either_nan)
                    correlation_matrix[i][j] = (same_measure_count / total_positions if total_positions > 0 else np.nan)
    
    correlation_df = pd.DataFrame(correlation_matrix, columns=marker_ratio_df.columns, index=marker_ratio_df.columns)
    correlation_df_percent = correlation_df * 100
    
    plt.figure(figsize=(25, 25))
    sns.heatmap(correlation_df_percent.iloc[::-1], annot=True, cmap='Reds', fmt=".0f", annot_kws={"size": 10})
    plt.title("%s SNP Concordance Heatmap" % runid, fontsize=20)
    plt.tight_layout()
    snp_output_heatmap = f'{rundir}/clinical/snpqc/{runid}_pct_new_markers_heatmap.png'
    plt.savefig(snp_output_heatmap)


# Main function with argparse
def main():
    parser = argparse.ArgumentParser(description="Process pileup files and generate marker ratio and heatmap.")
    parser.add_argument('--pileup_path', required=True, help='Path to the pileup files directory')
    parser.add_argument('--runid', required=True, help='Run ID for labeling the output')
    parser.add_argument('--rundir', required=True, help='Run output directory')
    parser.add_argument('--marker_bedfile', required=True, help='SNP marker bed file')
    args = parser.parse_args()
    
    marker_ratio_df = create_marker_ratiofile(args.pileup_path,args.rundir,args.marker_bedfile)
    marker_ratio_df_for_correlation = marker_ratio_df.drop(['POS', 'REF'], axis=1)
    create_correlation_matrix_heatmap(marker_ratio_df_for_correlation, args.runid, args.rundir)

if __name__ == "__main__":
    main()
