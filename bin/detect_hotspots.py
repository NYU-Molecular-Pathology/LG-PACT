#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import fnmatch
from typing import List, Tuple


def find(pattern, path):
    exclude_directories = set(['work'])
    file_path = []
    for root, dirs, files in os.walk(path):
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                file_path.append(os.path.join(root, name))
    return file_path


def output_caller_hotspots(out_df, caller, curr_sam, rundir):
    """Write the intermediate per-sample TSV file for Debugging"""
    fiName = curr_sam + "_Hotspots.tsv"
    outTsv = os.path.join(rundir, fiName)
    if os.path.isfile(outTsv):
        df1 = pd.read_csv(outTsv, dtype=str, sep='\t')
        df1[caller] = out_df[caller]
        df1.to_csv(outTsv, sep="\t", index=False)
    else:
        out_df.to_csv(outTsv, sep="\t", index=False)


def output_hotspots(out_df, pactid, rundir_output):
    """Write the final TSV output file"""
    fiName = pactid + "_Hotspots.tsv"
    outTsv = os.path.join(rundir_output, fiName)
    out_df.to_csv(outTsv, sep="\t", mode='a', index=False)
    print(f"File saved: {outTsv}")
    print("Hotspot report is rendered!")


def read_vcf(path: str) -> pd.Series:
    """Read a VCF file and return specific rows as a pd.Series."""
    try:
        df = pd.read_csv(path, comment='#', dtype=str, sep='\t', header=None, usecols=[0, 1])
        df = df[df[0].isin(["chr12", "chr7", "chr15", "chr3", "chr2"])]
        return (df[0] + ":" + df[1])
    except Exception as e:
        print("MuTect2:", e, "using columns=[0, 1] instead")
        df = pd.DataFrame(columns=[0, 1])
        return (df[0] + ":" + df[1])


def list_full_paths(directory: str) -> List[str]:
    """Return full paths of files in a directory."""
    return [entry.path for entry in os.scandir(directory) if entry.is_file()]


def get_sample_list(vcf_path: str) -> List[str]:
    """Extract sample list from vcf_path."""
    return list({f"{sam[6]}_{sam[7]}" for sam in (i.split('_') for i in vcf_path)})


def get_file_list(ts_number: List[str], vcf_path: List[str]) -> List[List[int]]:
    """Get indices of files matching ts_number in vcf_path."""
    return [[idx for idx, cont in enumerate(vcf_path) if substr in cont] for substr in ts_number]


def make_base_df(samList: List[str], temp_df: pd.DataFrame, hotspot_samples: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """Create base dataframe from samList and temp_df, dropping any extra samples if necessary."""

    # Ensure samList and hotspot_samples are paired correctly
    min_length = min(len(samList), len(hotspot_samples))
    dropped_sams = samList[min_length:] if len(samList) > min_length else []
    if dropped_sams:
        print(f"The following samples were dropped because they do not have corresponding NGS number: {dropped_sams}")
    paired_list = zip(samList[:min_length], hotspot_samples[:min_length])
    base_df = pd.concat([temp_df.assign(Sample=sam, Test_Number=hotspot) for sam, hotspot in paired_list], ignore_index=True)

    return base_df, dropped_sams


def sample_name(annotation_file: pd.DataFrame, column: str) -> List[str]:
    """Extract sample name from annotation_file's column."""
    return annotation_file[column].str.split('_', expand=True, n=5)[5].tolist()


def merge_vcf_col(caller_dict, hotspots_csv, out_df, ts_number, hotspot_samples):
    """_summary_
        Takes in template csv list of positions and checks each vcf file per sample for match
    Args:
        caller_dict (dict): dictionary of callers and their raw vcf file paths
        hotspots_csv (string): the path to the input csv file listing the fixed hotspot positions
        out_df (dataframe): the final output dataframe of YES and NO for detecting variants
        ts_number (list): list of tumor accession + '_' + dna_number
        hotspot_samples (list): list of the NGS-case numbers

    Returns:
        _type_: pandas.DataFrame
    """
    temp_df = pd.read_csv(hotspots_csv)
    out_df, dropped_samples = make_base_df(ts_number, temp_df, hotspot_samples)
    ts_number = [sample for sample in ts_number if sample not in dropped_samples]  # Exclude dropped samples

    pos_col = temp_df.columns.get_loc("Position")
    positions = temp_df.iloc[:, pos_col].tolist()

    for caller, inputDir in caller_dict.items():
        vcf_path = list(list_full_paths(inputDir))
        vcf_idx = get_file_list(ts_number, vcf_path)
        print(f"Checking caller: {caller}...")
        curr_col = out_df.columns.get_loc(caller)
        for idx, curr_sam in zip(vcf_idx, ts_number):
            print(f"Sample: {curr_sam}")
            filenames = [vcf_path[i] for i in idx if i < len(vcf_path)]  # Prevent index error
            if not filenames:  # If filenames list is empty, skip to the next sample
                print(f"No VCF files found for Sample: {curr_sam}")
                continue
            try:
                vcf_df_set = pd.concat(map(read_vcf, filenames)).tolist()
                yes_no = ['YES' if pos in vcf_df_set else 'NO' for pos in positions]
                curr_rows = out_df.index[out_df['Sample'] == curr_sam].tolist()
                out_df.iloc[curr_rows, curr_col] = yes_no
            except ValueError as e:
                print(f"Error processing VCF files for Sample: {curr_sam}. Error: {e}")

    return out_df


def startParse(hotspot_samples, rundir_output, ts_number, rundir):
    out_df = None
    caller_dict = {
        "LoFreqSomatic": os.path.join(rundir_output, "variants/LoFreqSomatic/norm/"),
        "Strelka": os.path.join(rundir_output, "variants/Strelka/raw/"),
        "Mutect": os.path.join(rundir_output, "variants/MuTect2/raw/")
    }
    hotspot_find_file = find('hotspot_genes.csv', rundir)
    print(hotspot_find_file)
    hotspots_csv = hotspot_find_file[0]
    return merge_vcf_col(caller_dict, hotspots_csv, out_df, ts_number, hotspot_samples)


def main(runid, pactid, rundir):
    print("Starting detect_hotspots.py for " + runid)
    rundir_output = os.path.join(rundir, 'output/')
    demux_ss_find = rundir + 'demux-samplesheet.csv'
    print("demux_ss_find:" + demux_ss_find)
    demux_ss_file = pd.read_csv(demux_ss_find, skiprows=19)
    is_tumor = demux_ss_file['Tumor_Content'] != 0
    ts_number = demux_ss_file.loc[is_tumor, 'Specimen_ID'].tolist()
    hotspot_samples = demux_ss_file['Test_Number'].unique()
    hotspot_samples = hotspot_samples[hotspot_samples != '0']
    hotspot_variants = startParse(hotspot_samples, rundir_output, ts_number, rundir)
    output_hotspots(hotspot_variants, pactid, rundir_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory's output dir")
    parser.add_argument('-rid', '--runid', help='Your run id', required=True)
    parser.add_argument('-pactid', '--pactid', help='PACT ID', required=True)
    parser.add_argument('-rdir', '--rundir', help='Your run dir', required=True)
    args = parser.parse_args()
    main(args.runid, args.pactid, args.rundir)