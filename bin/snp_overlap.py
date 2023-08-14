#!/usr/bin/env python

__author__ = "Kelsey Zhu"
__version__ = "1.0.1"

import pandas as pd
import os
import argparse

def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="NGS607 output directory")
    parser.add_argument("-p", "--percent", type=float, required=False,
                        help="snp overlaping percentage",
                        default=0.80)
    return parser.parse_args()

def get_snp_overlap(tsv_path, tumor, normal, cutoff):
    tumor_snps = pd.read_csv(os.path.join(tsv_path, tumor), sep='\t', skiprows=0)
    tumor_snps = tumor_snps.loc[tumor_snps['FREQ'].astype(float) >= cutoff]
    normal_snps = pd.read_csv(os.path.join(tsv_path, normal), sep='\t', skiprows=0)
    normal_snps = normal_snps.loc[normal_snps['FREQ'].astype(float) >= cutoff]
    tumor_snps['varkey'] = tumor_snps['CHROM'] + ":" + tumor_snps['POS'].astype(str) + ":" + tumor_snps['REF'] + ">" + \
                           tumor_snps['ALT']
    print("total tumor snps: %s" %len(tumor_snps.index))
    normal_snps['varkey'] = normal_snps['CHROM'] + ":" + normal_snps['POS'].astype(str) + ":" + normal_snps[
        'REF'] + ">" + normal_snps['ALT']
    print("total normal snps (denominator): %s" %len(normal_snps.index))
    common_snps = tumor_snps.merge(normal_snps, left_on='varkey', right_on='varkey', how='inner')
    print("overlap snps: %s" %len(common_snps.index))
    combined_snps = pd.concat([tumor_snps, normal_snps], axis=0)
    combined_snps = combined_snps.drop_duplicates(subset='varkey', keep="first")
    print("total unique snps: %s" %len(combined_snps.index))
    if len(normal_snps.index) != 0:
        snp_pcent = round((len(common_snps.index) / len(normal_snps.index))*100,2)
    else:
        snp_pcent = 0
    print("snp overlap: %s" % snp_pcent)
    return str(len(tumor_snps.index)), str(len(normal_snps.index)), str(len(common_snps.index)), str(snp_pcent)

def main():
    args = get_options()
    work_dir = args.output
    percent = args.percent
    tsv_path = os.path.join(work_dir, "variants/LoFreq/tsv")
    with open(os.path.join(work_dir,'samplesheet.tsv')) as f, \
            open(os.path.join(work_dir, "snp_overlap_%s.tsv"%str(int(percent*100))),'w') as of:
        of.write("Tumor\tNormal\tT_SNP_Count\tN_SNP_Count\tSNP_Overlap_Count\tSNP_Overlap_Pcent(%)\n")
        for line in f.readlines()[1:]:
            parts = line.split("\t")
            if parts[2] == "NA": continue
            t_snp_count, n_snp_count, overlap_snp, snp_overlap_pcent = get_snp_overlap(tsv_path,
                        "%s.LoFreq.NA.reformat.tsv"%parts[1], "%s.LoFreq.NA.reformat.tsv"%parts[2], percent)
            of.write("%s\n"%"\t".join([parts[1],parts[2],t_snp_count,n_snp_count,overlap_snp, snp_overlap_pcent]))

if __name__ == "__main__":
    main()

