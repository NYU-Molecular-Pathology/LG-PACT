#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calcalate the tumor mutation burden in variants/Megabase for Validation

"""
import os
import csv
import argparse
import pandas as pd

"""
Additional criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5%
2- coverage > 200X.
3- include Non Synonymous and synonymous SNV.
4- include only exonic.

For Unmatched. Add condition below;
5- exclude Cosmic HS
6- Frequency of <0.4% in ExAC.
"""

NA_strs = ['.'] # strings used as "NA" values in the table
Func_refGene_allowed = ['exonic'] # , 'splicing' 'exonic;splicing', , 'UTR5'
coverage_min = 200.0 # should correspond to GATK CallableLoci depth cutoff
frequency_min = [0.05, 0.10]
ExAC_max = 0.004
ExonicFunc_refGene_allowed = ['nonsynonymous SNV','synonymous SNV']

def filter_rules(row, type="paired"):
    """
    Return True or False if the row passes all the filter criteria
    """
    frequency = float(row['AF'])
    coverage = float(row['DP'])
    Func_refGene = row['Func.refGene']
    ExonicFunc_refGene = row['ExonicFunc.refGene']
    if type != "paired":
        ExAC_value = row['ExAC_ALL']
    frequency_pass = frequency > frequency_min[0] if type == "paired" else frequency_min[1]
    coverage_pass = coverage >= coverage_min
    not_in_COSMIC = True
    if type != "paired":
        not_in_COSMIC = row['cosmic70'] in NA_strs
    in_Func_refGene_allowed = Func_refGene in Func_refGene_allowed
    in_Func_ExonicFunc_refGene = ExonicFunc_refGene in ExonicFunc_refGene_allowed
    in_ExAC_allowed = True
    if type == "unpaired":
        in_ExAC_allowed = ExAC_value in NA_strs or float(ExAC_value) <= ExAC_max

    return(all([in_Func_refGene_allowed, not_in_COSMIC, coverage_pass,
                frequency_pass, in_ExAC_allowed, in_Func_ExonicFunc_refGene]))

def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--loci", type=str, required=True,
                        help="callable loci from alignments")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to input variant annotation including file name")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output TMB file including file name")
    parser.add_argument("-s", "--sample_sheet", type=str, required=True, help="sample sheet")
    parser.add_argument("-t", "--type", type=str, required=False,
                        help="tumor sample type, paired not unpaired", default="paired")
    return parser.parse_args()

def main(callable_loci_file,annotation_paired_file,demux_ss,callerType,output_tmb_file):
    loci_dict = dict()
    tmb_dict = dict()
    sample_dict = dict()
    excl_list = ["SC-SERACARE","NC-HAPMAP","SC_SERACARE","NC_HAPMAP"]
    with open (callable_loci_file) as fl:
        for line in fl:
            items = line.strip().split()
            loci_dict[items[1]] = items[0]
    with open(demux_ss) as fs:
        for line in fs.readlines()[20:]:
            items = line.split(",")
            if items[0].endswith(tuple(excl_list)) or items[2] == "NA": continue
            sample_dict[items[0]] = items[-7].strip()
            print(items[0], items[-7].strip())
    with open(annotation_paired_file) as fin, open(output_tmb_file, 'w') as fout:
        fout.write("SampleID\tVariantCaller\tnBases\tnVariants\tTMB\tTumorType\n")
        reader = csv.DictReader(fin, delimiter = '\t')
        print(reader.fieldnames)
        for row in reader:
            if filter_rules(row):
                caller = row['VariantCaller']
                if caller not in tmb_dict.keys():
                    tmb_dict[caller] = {}
                sample = row['Tumor'] if callerType == "paired" else row['Sample']
                if sample.endswith(tuple(excl_list)) : continue
                if sample not in tmb_dict[caller].keys():
                    tmb_dict[caller][sample] = {"variants":0}
                tmb_dict[caller][sample]['variants'] += 1

        #print (tmb_dict)
        for caller in tmb_dict.keys():
            for sample in tmb_dict[caller].keys():
                val = tmb_dict[caller][sample]
                tmb = round(float(val['variants'])/float(loci_dict[sample])*1000000,2)
                fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,caller,str(loci_dict[sample]),
                                                       str(val['variants']),str(tmb), sample_dict[sample]))

if __name__ == "__main__":
    main()