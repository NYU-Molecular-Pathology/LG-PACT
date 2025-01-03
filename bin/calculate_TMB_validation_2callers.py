#!/usr/bin/env python3
"""
Calculate the tumor mutation burden in variants/Megabase for Validation (MuTect2+LoFreqSomatic)

"""
import argparse
import pandas as pd
import os

"""
Additional criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5%
2- coverage > 200X.
3- include Non Synonymous, synonymous SNV, stoploss and stopgain
4- include only exonic.
5- variants called by both MuTect and LoFreqSomatic only

"""
def load_data(demux_ss, annot_paired):
    """
    Load and parse the demultiplexing and annotation data.

    Parameters:
    - demux_ss (str): Path to the CSV file containing demultiplexing data.
    - annot_paired (str): Path to the TSV file containing paired annotation data.

    Returns:
    - tuple: A tuple containing:
        - demux_ss_df (pd.DataFrame): DataFrame with demultiplexing data.
        - annot_paired_df (pd.DataFrame): DataFrame with annotation data.
    """
    # Load demultiplexing data, skipping the first 19 rows for metadata
    demux_ss_df = pd.read_csv(demux_ss, skiprows=19)
    # Load annotation data with tab-separated format
    annot_paired_df = pd.read_csv(annot_paired, sep="\t")
    return demux_ss_df, annot_paired_df

def filter_annotation_data(annot_paired_df):
    """
    Filter annotation data based on specific criteria.

    Parameters:
    - annot_paired_df (pd.DataFrame): DataFrame containing annotation data.

    Returns:
    - pd.DataFrame: Filtered annotation data including only relevant variants.
    """
    # Convert the 'POS' column to string type
    annot_paired_data = annot_paired_df.astype({'POS': str})
    # Filter rows with AF > 0.05 and DP >= 200
    annot_paired_filter_AFDP = annot_paired_data[(annot_paired_data['AF'] > 0.05) & (annot_paired_data['DP'] >= 200)]
    # Further filter for exonic variants with specific functional annotations
    annot_paired_type_filter = annot_paired_filter_AFDP[((annot_paired_filter_AFDP['Func.refGene'] == "exonic") & (annot_paired_filter_AFDP['ExonicFunc.refGene'].isin(["nonsynonymous SNV", "synonymous SNV", "stoploss", "stopgain"])))]
    # Create a new 'Variant' column by concatenating CHROM, POS, REF, and ALT columns
    annot_paired_type_filter['Variant'] = annot_paired_type_filter[['CHROM', 'POS', 'REF', 'ALT']].apply(lambda x: ':'.join(x), axis=1)
    # Filter variants called by MuTect2 or LoFreqSomatic
    annot_Mutect_LoFreq = annot_paired_type_filter[(annot_paired_type_filter['VariantCaller'].isin(["MuTect2", "LoFreqSomatic"]))]
    # Exclude control samples based on their Tumor column values
    annot_Mutect_LoFreq_nocontrols = annot_Mutect_LoFreq[~annot_Mutect_LoFreq['Tumor'].str.contains("SC-SERACARE|NC-HAPMAP|SC_SERACARE|NC_HAPMAP")]
    return annot_Mutect_LoFreq_nocontrols

def filter_demux_data(demux_ss_df):
    """
    Filter and process demultiplexing data.

    Parameters:
    - demux_ss_df (pd.DataFrame): DataFrame containing demultiplexing data.

    Returns:
    - tuple: A tuple containing:
        - demux_ss_data_final (pd.DataFrame): Filtered demultiplexing data.
        - demux_ss_dict (dict): Dictionary mapping sample names to tumor types.
    """
    # Drop rows with missing values
    demux_ss_data = demux_ss_df.dropna()
    # Exclude control samples based on their Sample_Name column values
    demux_ss_data_final = demux_ss_data[~demux_ss_data['Sample_Name'].str.contains("SC_SERACARE|NC_HAPMAP|NTC_H20")]
    # Create a dictionary mapping Sample_Name to Tumor_Type
    demux_ss_dict = dict(zip(demux_ss_data_final.Sample_Name, demux_ss_data_final.Tumor_Type))
    return demux_ss_data_final, demux_ss_dict

def compute_variant_counts(samples, annot_df, output_dir):
    """
    Computes common variant counts for each sample and writes annotated variants to Excel files.

    Parameters:
    - samples (list): List of sample identifiers.
    - annot_df (pd.DataFrame): DataFrame containing annotated variants.
    - output_dir (str): Directory to save output Excel files.

    Returns:
    - dict: Dictionary mapping each sample to the number of common variants.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Columns to retain in the output
    selected_cols = ["CHROM", "POS", "REF", "ALT", "DP", "AF", "VariantCaller",
                     "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"]
    # Initialize the result dictionary
    sample_nvariants_dict = {}
    for sample in samples:
        # Filter DataFrame for the current sample
        sample_df = annot_df[annot_df['Tumor'] == sample]
        # Separate MuTect2 and LoFreqSomatic variants
        mutect_variants = sample_df[sample_df['VariantCaller'] == "MuTect2"]
        lofreq_variants = sample_df[sample_df['VariantCaller'] == "LoFreqSomatic"]
        # Find common variants
        common_variants = pd.merge(mutect_variants, lofreq_variants, how='inner', on='Variant', suffixes=('_mutect', '_lfs'))
        # Prepare DataFrames for writing
        mutect2_df = common_variants.filter(regex='_mutect').rename(columns=lambda x: x.replace("_mutect", ""))
        lofreqsomatic_df = common_variants.filter(regex='_lfs').rename(columns=lambda x: x.replace("_lfs", ""))
        # Write to Excel
        output_path = os.path.join(output_dir, f"{sample}.xlsx")
        with pd.ExcelWriter(output_path) as writer:
            mutect2_df[selected_cols].to_excel(writer, sheet_name="MuTect2", index=False)
            lofreqsomatic_df[selected_cols].to_excel(writer, sheet_name="LoFreqSomatic", index=False)
        # Update the dictionary with the count of common variants
        sample_nvariants_dict[sample] = len(common_variants)
    return sample_nvariants_dict

def load_callable_loci(callable_loci):
    """
    Load callable loci data into a dictionary.

    Parameters:
    - callable_loci (str): Path to the callable loci file.

    Returns:
    - dict: Dictionary mapping loci to their corresponding values.
    """
    loci_dict = {}
    with open(callable_loci) as fl:
        for line in fl:
            # Split each line and map loci (second column) to its value (first column)
            items = line.strip().split()
            loci_dict[items[1]] = int(items[0])
    return loci_dict

def merge_dicts(sample_nvariants_dict, loci_dict, demux_ss_dict):
    """
    Merge multiple dictionaries to create a combined data structure.

    Parameters:
    - sample_nvariants_dict (dict): Dictionary of sample IDs to variant counts.
    - loci_dict (dict): Dictionary of loci and their associated values.
    - demux_ss_dict (dict): Dictionary of sample IDs to tumor types.

    Returns:
    - dict: Dictionary mapping sample IDs to a tuple of merged values.
    """
    # Combine dictionaries into a single dictionary with tuples of values
    dict_list = [sample_nvariants_dict, loci_dict, demux_ss_dict]
    tmb_vals_dict = {k: tuple(d[k] for d in dict_list) for k in sample_nvariants_dict.keys()}
    return tmb_vals_dict

def write_output(tmb_vals_dict, output_tmb_file):
    """
    Write TMB values to an output file.

    Parameters:
    - tmb_vals_dict (dict): Dictionary of sample IDs to their calculated TMB values.
    - output_tmb_file (str): Path to the output file.

    Returns:
    - None
    """
    # Open the output file for writing
    with open(output_tmb_file, 'w') as fout:
        # Write header line
        fout.write("SampleID\tVariantCaller\tnBases\tnVariants\tTMB\tTumorType\n")
    
        for sampleID, values in tmb_vals_dict.items():
            
            nBases = values[1]
            nVariants = values[0]
            tumorType = values[2]

            # Determine TMB and nVariants outputs
            tmb = "NA" if nBases == 0 or nBases < 1000000 else round(nVariants / nBases * 1000000, 2)
            nVariants_output = "NA" if nBases < 1000000 else nVariants
            # Write the line
            fout.write(f"{sampleID}\tMuTect2+LoFreqsomatic\t{nBases}\t{nVariants_output}\t{tmb}\t{tumorType}\n")
               
def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--loci", type=str, required=True,
                        help="callable loci from alignments")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to input variant annotation including file name")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output TMB file including file name")
    parser.add_argument("-s", "--sample_sheet", type=str, required=True, help="sample sheet")
    parser.add_argument("-vao","--variant_annotation_output",type=str, required=True, help="output path for TMB variants")
    parser.add_argument("-t", "--type", type=str, required=False,
                        help="tumor sample type, paired not unpaired", default="paired")
    return parser.parse_args()

def main():
    args = get_options()
    # Load demux ss and annotations paired data
    demux_ss_df, annot_paired_df = load_data(args.sample_sheet, args.input)
    # Apply filter criteria for TMB
    annot_Mutect_LoFreq_nocontrols = filter_annotation_data(annot_paired_df)
    # Get sample and tumor type list
    demux_ss_data_final, demux_ss_dict = filter_demux_data(demux_ss_df)
    samples = demux_ss_data_final.Sample_Name.unique()
    # Get variant count for each sample after applying filters
    sample_nvariants_dict = compute_variant_counts(samples, annot_Mutect_LoFreq_nocontrols, args.variant_annotation_output)
    # Load the callable loci data
    loci_dict = load_callable_loci(args.loci)
    # Calculate TMB
    tmb_vals_dict = merge_dicts(sample_nvariants_dict, loci_dict, demux_ss_dict)
    # Write final TMB output
    write_output(tmb_vals_dict,args.output)

if __name__ == "__main__":
    main()