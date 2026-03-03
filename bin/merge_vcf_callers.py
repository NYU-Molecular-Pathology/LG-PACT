#!/usr/bin/env python3
"""
Merge VCF files from MuTect, LoFreqSomatic, and Strelka somatic variant callers.

Priority rules for selecting which caller's INFO/FORMAT values to use:
1. If variant present in all 3 callers -> use MuTect INFO/FORMAT
2. If variant in MuTect + Strelka only -> use MuTect INFO/FORMAT
3. If variant in LoFreqSomatic + Strelka only -> use LoFreqSomatic INFO/FORMAT
4. Otherwise -> use the caller that detected it

Additional logic:
- Within each priority rule, preferentially use values from callers with FILTER=PASS
- SOURCE tag shows which callers had FILTER=PASS for each variant
- For non-PASS variants, SOURCE shows all callers that detected it
"""

import sys
import argparse
import gzip
import os
from collections import defaultdict
from pathlib import Path

__version__ = "2.0.0"


def open_vcf(filename):
    """Open VCF file (handles gzipped files)."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')


def read_vcf(filename):
    """Read VCF file and return header lines and variants."""
    header_lines = []
    variants = {}
    
    with open_vcf(filename) as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line.rstrip())
            elif line.startswith('#CHROM'):
                header_lines.append(line.rstrip())
            else:
                fields = line.rstrip().split('\t')
                chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
                var_key = f"{chrom}:{pos}:{ref}:{alt}"
                variants[var_key] = fields
    
    return header_lines, variants


def get_sample_count(header):
    """Get number of sample columns from header."""
    for line in header:
        if line.startswith('#CHROM'):
            fields = line.split('\t')
            if len(fields) > 9:
                return len(fields) - 9
    return 0


def normalize_variant_line(fields, expected_samples):
    """Ensure variant has correct number of sample columns."""
    if len(fields) <= 8 and expected_samples > 0:
        fields.append('GT')
        for _ in range(expected_samples):
            fields.append('./.')
    elif len(fields) == 9:
        for _ in range(expected_samples):
            fields.append('./.')
    elif len(fields) < 9 + expected_samples:
        while len(fields) < 9 + expected_samples:
            fields.append('./.')
    return fields


def merge_headers(mutect_header, lofreqsomatic_header, strelka_header):
    """Merge headers from all three VCF files."""
    merged = []
    seen_lines = set()
    
    for line in mutect_header:
        if line.startswith('##fileformat') or line.startswith('##reference'):
            merged.append(line)
            seen_lines.add(line)
    
    merged.append(f'##source=VCFMerger_v{__version__}')
    
    meta_lines = []
    for header in [mutect_header, lofreqsomatic_header, strelka_header]:
        for line in header:
            if any(line.startswith(f'##{tag}=') for tag in ['INFO', 'FORMAT', 'FILTER', 'contig', 'ALT']):
                if line not in seen_lines:
                    meta_lines.append(line)
                    seen_lines.add(line)
    
    merged.extend(sorted(meta_lines))
    
    for line in mutect_header:
        if line.startswith('#CHROM'):
            merged.append(line)
            break
    
    return merged


def merge_variants(mutect_vars, lofreqsomatic_vars, strelka_vars):
    """Merge variants with FILTER-based prioritization."""
    merged = {}
    all_keys = set(mutect_vars.keys()) | set(lofreqsomatic_vars.keys()) | set(strelka_vars.keys())
    
    for var_key in all_keys:
        in_mutect = var_key in mutect_vars
        in_lofreqsomatic = var_key in lofreqsomatic_vars
        in_strelka = var_key in strelka_vars
        
        mutect_is_pass = in_mutect and mutect_vars[var_key][6] == 'PASS'
        lofreq_is_pass = in_lofreqsomatic and lofreqsomatic_vars[var_key][6] == 'PASS'
        strelka_is_pass = in_strelka and strelka_vars[var_key][6] == 'PASS'
        
        # Build SOURCE tag with PASS callers
        source_parts = []
        if mutect_is_pass:
            source_parts.append('Mutect')
        if lofreq_is_pass:
            source_parts.append('LoFreqSomatic')
        if strelka_is_pass:
            source_parts.append('Strelka')
        
        # If no PASS callers, use all callers that detected it
        if source_parts:
            source = '+'.join(source_parts)
        else:
            nonpass_parts = []
            if in_mutect:
                nonpass_parts.append('Mutect')
            if in_lofreqsomatic:
                nonpass_parts.append('LoFreqSomatic')
            if in_strelka:
                nonpass_parts.append('Strelka')
            source = '+'.join(nonpass_parts)
        
        # Select which caller's values to use (priority: MuTect > LoFreq > Strelka)
        if mutect_is_pass:
            merged[var_key] = (mutect_vars[var_key], source)
        elif lofreq_is_pass:
            merged[var_key] = (lofreqsomatic_vars[var_key], source)
        elif strelka_is_pass:
            merged[var_key] = (strelka_vars[var_key], source)
        elif in_mutect:
            merged[var_key] = (mutect_vars[var_key], source)
        elif in_lofreqsomatic:
            merged[var_key] = (lofreqsomatic_vars[var_key], source)
        elif in_strelka:
            merged[var_key] = (strelka_vars[var_key], source)
    
    return merged


def sort_variants(merged_vars):
    """Sort variants by chromosome and position."""
    def get_sort_key(item):
        var_key, (fields, source) = item
        chrom = fields[0]
        pos = int(fields[1])
        
        chrom_clean = chrom.replace('chr', '')
        if chrom_clean.isdigit():
            chrom_sort = int(chrom_clean)
        elif chrom_clean == 'X':
            chrom_sort = 23
        elif chrom_clean == 'Y':
            chrom_sort = 24
        elif chrom_clean in ['M', 'MT']:
            chrom_sort = 25
        else:
            chrom_sort = 26
        
        return (chrom_sort, pos)
    
    return sorted(merged_vars.items(), key=get_sort_key)


def main():
    parser = argparse.ArgumentParser(
        description='Merge VCF files from MuTect2, LoFreqSomatic, and Strelka',
        epilog='Example: python merge_vcf.py -m mutect.vcf -l lofreq.vcf -s strelka.vcf -o merged.vcf --add-source-tag'
    )
    
    parser.add_argument('-m', '--mutect', required=True, help='MuTect2 VCF file')
    parser.add_argument('-l', '--lofreqsomatic', required=True, help='LoFreqSomatic VCF file')
    parser.add_argument('-s', '--strelka', required=True, help='Strelka VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output merged VCF file')
    parser.add_argument('--add-source-tag', action='store_true',
                        help='Add SOURCE tag showing which callers had FILTER=PASS')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    print("Reading VCF files...")
    mutect_header, mutect_vars = read_vcf(args.mutect)
    lofreqsomatic_header, lofreqsomatic_vars = read_vcf(args.lofreqsomatic)
    strelka_header, strelka_vars = read_vcf(args.strelka)
    
    print(f"  MuTect: {len(mutect_vars)} variants")
    print(f"  LoFreq: {len(lofreqsomatic_vars)} variants")
    print(f"  Strelka: {len(strelka_vars)} variants")
    
    print("Merging variants...")
    merged_header = merge_headers(mutect_header, lofreqsomatic_header, strelka_header)
    expected_samples = get_sample_count(merged_header)
    
    if args.add_source_tag:
        source_line = '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Variant caller(s) where FILTER=PASS">'
        merged_header.insert(-1, source_line)
    
    merged_vars = merge_variants(mutect_vars, lofreqsomatic_vars, strelka_vars)
    sorted_vars = sort_variants(merged_vars)
    
    print(f"Writing {len(sorted_vars)} variants to {args.output}...")
    with open(args.output, 'w') as out:
        for line in merged_header:
            out.write(line + '\n')
        
        for var_key, (fields, source) in sorted_vars:
            fields = normalize_variant_line(fields, expected_samples)
            
            if args.add_source_tag:
                info_field = fields[7]
                if info_field == '.':
                    fields[7] = f'SOURCE={source}'
                else:
                    fields[7] = f'{info_field};SOURCE={source}'
            
            out.write('\t'.join(fields) + '\n')
    
    print("\nSummary:")
    sources = defaultdict(int)
    for _, (_, source) in sorted_vars:
        sources[source] += 1
    
    for source in sorted(sources.keys()):
        print(f"  {source}: {sources[source]} variants")
    
    print(f"\nDone! Output: {args.output}")


if __name__ == '__main__':
    main()