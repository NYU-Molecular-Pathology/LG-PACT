#!/usr/bin/env python3

import os
import pandas as pd
import jinja2
import argparse
import fnmatch
from typing import List
import pandas as pd


TEMPLATE_FILE_HOTSPOT = "hotspot_template.html"

def find(pattern, path):
    exclude_directories = set(['work'])
    file_path = []
    for root,dirs,files in os.walk(path):
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


def output_hotspots(out_df,  pactid, rundir_output):
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
    except:
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


def make_base_df(samList: List[str], temp_df: pd.DataFrame, hotspot_samples: List[str]) -> pd.DataFrame:
    """Create base dataframe from samList and temp_df."""
    return pd.concat([temp_df.assign(Sample=curr_sam, Test_Number=hotspot_samples[idx]) for idx, curr_sam in enumerate(samList)], ignore_index=True)


def sample_name(annotation_file: pd.DataFrame, column: str) -> List[str]:
    """Extract sample name from annotation_file's column."""
    return annotation_file[column].str.split('_', expand=True, n=5)[5].tolist()


def merge_vcf_col(caller_dict, hotspots_csv, out_df, ts_number, hotspot_samples, rundir_output):
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
    out_df = make_base_df(ts_number, temp_df, hotspot_samples)
    pos_col = temp_df.columns.get_loc("Position")
  
    positions = temp_df.iloc[:, pos_col].tolist()
#    print("positions", positions)
    
    for caller, inputDir in caller_dict.items():
        vcf_path = list(list_full_paths(inputDir))
        vcf_idx = get_file_list(ts_number, vcf_path)
        print(f"Checking caller: {caller}...")
        curr_col = out_df.columns.get_loc(caller)
        for idx, curr_sam in zip(vcf_idx, ts_number):
            print(f"Sample: {curr_sam}")
            filenames = [vcf_path[i] for i in idx]
            vcf_df_set = pd.concat(map(read_vcf, filenames)).tolist()
 #           print("vcf_df_set\n", vcf_df_set)
            yes_no = ['YES' if pos in vcf_df_set else 'NO' for pos in positions]
            curr_rows = out_df.index[out_df['Sample'] == curr_sam].tolist()
  #          print("curr_rows", curr_rows)
            out_df.iloc[curr_rows, curr_col] = yes_no
            
    return out_df
   
def startParse(hotspot_samples, rundir_output, ts_number, runid, rundir):
    out_df = None
    caller_dict = {
        "LoFreqSomatic": os.path.join(rundir_output, "variants/LoFreqSomatic/norm/"),
        "Strelka": os.path.join(rundir_output, "variants/Strelka/raw/"),
        "Mutect": os.path.join(rundir_output, "variants/MuTect2/raw/")
    }
    maindir = rundir
    hotspot_find_file = find('hotspot_genes.csv', maindir)
    hotspots_csv = hotspot_find_file[0]
    return merge_vcf_col(caller_dict, hotspots_csv, out_df, ts_number, hotspot_samples, rundir_output)

def make_html_report(runid: str, rundir: str, pactid: str, hotspot_variants: pd.DataFrame, hotspot_samples: List[str]):
    template_path = os.path.join(rundir,"template/")
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    html_out = os.path.join(rundir, "clinical")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template_hotspot = templateEnv.get_template(TEMPLATE_FILE_HOTSPOT)
    html_out_hotspot = template_hotspot.render(hotspot_variants=hotspot_variants, hotspot_samples=hotspot_samples, pactid=pactid, runid=runid)
    os.makedirs(html_out, exist_ok=True)
    hotspot_filename = os.path.join(html_out, f'{pactid}-hotspot_Variants.html')
    with open(hotspot_filename, 'w') as html_file_hotspot:
        html_file_hotspot.write(html_out_hotspot)
    os.system(f"chmod -R g+rwx {html_out}")
    print("Hotspot report is rendered!")

def main(runid, pactid, rundir):
    rundir_output = os.path.join(rundir, 'output/')
    maindir = rundir
    demux_ss_find = find('demux-samplesheet.csv',maindir)
    demux_ss_file = pd.read_csv(demux_ss_find[0], skiprows=19)
    is_tumor = demux_ss_file['Tumor_Content'] != 0
    ts_number = demux_ss_file.loc[is_tumor, 'Specimen_ID'].tolist()
    hotspot_samples = demux_ss_file['Test_Number'].unique()
    hotspot_samples = hotspot_samples[hotspot_samples != '0']
    hotspot_variants = startParse(hotspot_samples, rundir_output, ts_number, runid, rundir)
    output_hotspots(hotspot_variants, pactid, rundir_output)
    #make_html_report(runid, rundir, pactid, hotspot_variants, hotspot_samples)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory's output dir")
    parser.add_argument('-rid', '--runid', help='Your run id', required=True)
    parser.add_argument('-pactid', '--pactid', help='PACT ID', required=True)
    parser.add_argument('-rdir','--rundir', help='Your run dir', required=True )
    args = parser.parse_args()
    main(args.runid, args.pactid, args.rundir)

