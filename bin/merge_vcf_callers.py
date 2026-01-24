#!/usr/bin/env python3
"""
Merge VCF files from MuTect, LoFreqSomatic, and Strelka somatic variant callers.

Priority rules:
1. If variant present in all 3 callers -> use MuTect INFO/FORMAT
2. If variant in MuTect + Strelka only -> use MuTect INFO/FORMAT
3. If variant in LoFreqSomatic + Strelka only -> use LoFreqSomatic INFO/FORMAT
4. Otherwise -> use the caller that detected it
"""

import sys
import argparse
from collections import defaultdict
import gzip


def open_vcf(filename):
    """Open VCF file (handles gzipped files)."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')


def parse_vcf_line(line):
    """Parse a VCF data line into components."""
    fields = line.strip().split('\t')
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    
    # Create variant key (chrom:pos:ref:alt)
    var_key = f"{chrom}:{pos}:{ref}:{alt}"
    
    return var_key, fields


def read_vcf(filename, caller_name):
    """
    Read VCF file and return header lines and variant dictionary.
    
    Returns:
        header_lines: list of header lines (including #CHROM line)
        variants: dict mapping var_key to full VCF line fields
    """
    header_lines = []
    variants = {}
    
    with open_vcf(filename) as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line.rstrip())
            elif line.startswith('#CHROM'):
                header_lines.append(line.rstrip())
            else:
                var_key, fields = parse_vcf_line(line)
                variants[var_key] = fields
    
    print(f"Read {len(variants)} variants from {caller_name}", file=sys.stderr)
    return header_lines, variants


def get_sample_count(header):
    """Get the number of sample columns from the header."""
    for line in header:
        if line.startswith('#CHROM'):
            fields = line.split('\t')
            # Fields after FORMAT are sample columns
            if len(fields) > 9:
                return len(fields) - 9
            else:
                return 0
    return 0


def normalize_variant_line(fields, expected_samples):
    """
    Ensure variant line has the correct number of sample columns.
    Adds placeholder columns if needed.
    
    Args:
        fields: List of VCF fields for a variant
        expected_samples: Number of sample columns expected
    
    Returns:
        Normalized list of fields
    """
    # Standard VCF has 8 mandatory columns before FORMAT
    # If there's a FORMAT column, then sample columns follow
    current_fields = len(fields)
    
    # If we have fewer than 9 fields, no FORMAT/samples present
    if current_fields <= 8:
        if expected_samples > 0:
            # Add FORMAT column with minimal genotype
            fields.append('GT')
            # Add placeholder sample columns (missing genotype)
            for _ in range(expected_samples):
                fields.append('./.')
    elif current_fields == 9:
        # Has FORMAT but no samples - add placeholder samples
        for _ in range(expected_samples):
            fields.append('./.')
    elif current_fields < 9 + expected_samples:
        # Has some samples but not enough - add more placeholders
        while len(fields) < 9 + expected_samples:
            fields.append('./.')
    
    return fields


def merge_headers(mutect_header, lofreqsomatic_header, strelka_header):
    """
    Merge headers from all three VCF files.
    Keep all unique ##INFO, ##FORMAT, ##FILTER lines.
    Use MuTect header as base.
    """
    merged = []
    seen_lines = set()
    
    # Add fileformat and reference from MuTect
    for line in mutect_header:
        if line.startswith('##fileformat') or line.startswith('##reference'):
            merged.append(line)
            seen_lines.add(line)
            
    # Add source information
    merged.append('##source=MergedVCF_MuTect_LoFreqSomatic_Strelka')
    
    # Collect all INFO, FORMAT, FILTER, contig lines from all callers
    meta_lines = []
    for header in [mutect_header, lofreqsomatic_header, strelka_header]:
        for line in header:
            if any(line.startswith(f'##{tag}=') for tag in ['INFO', 'FORMAT', 'FILTER', 'contig', 'ALT']):
                if line not in seen_lines:
                    meta_lines.append(line)
                    seen_lines.add(line)
    
    merged.extend(sorted(meta_lines))
    
    # Add the column header line (from MuTect)
    for line in mutect_header:
        if line.startswith('#CHROM'):
            merged.append(line)
            break
    
    return merged


def merge_variants(mutect_vars, lofreqsomatic_vars, strelka_vars):
    """
    Merge variants according to priority rules.
    
    Returns:
        dict mapping var_key to (fields, source) where source is the caller used
    """
    merged = {}
    
    # Get all unique variant keys
    all_keys = set(mutect_vars.keys()) | set(lofreqsomatic_vars.keys()) | set(strelka_vars.keys())
    
    for var_key in all_keys:
        in_mutect = var_key in mutect_vars
        in_lofreqsomatic = var_key in lofreqsomatic_vars
        in_strelka = var_key in strelka_vars
        
        # Rule 1: If in all 3, use MuTect
        if in_mutect and in_lofreqsomatic and in_strelka:
            merged[var_key] = (mutect_vars[var_key], 'MuTect+LoFreqSomatic+Strelka')
        
        # Rule 2: If in MuTect + Strelka, use MuTect
        elif in_mutect and in_strelka and not in_lofreqsomatic:
            merged[var_key] = (mutect_vars[var_key], 'MuTect+Strelka')
        
        # Rule 3: If in LoFreqSomatic + Strelka only, use LoFreqSomatic
        elif in_lofreqsomatic and in_strelka and not in_mutect:
            merged[var_key] = (lofreqsomatic_vars[var_key], 'LoFreqSomatic+Strelka')
        
        # Rule 4: Use whichever caller detected it
        elif in_mutect:
            merged[var_key] = (mutect_vars[var_key], 'MuTect')
        elif in_lofreqsomatic:
            merged[var_key] = (lofreqsomatic_vars[var_key], 'LoFreqSomatic')
        elif in_strelka:
            merged[var_key] = (strelka_vars[var_key], 'Strelka')
    
    return merged


def sort_variants(merged_vars):
    """
    Sort variants by chromosome and position.
    """
    # Extract chrom and pos for sorting
    def get_sort_key(item):
        var_key, (fields, source) = item
        chrom = fields[0]
        pos = int(fields[1])
        
        # Handle chromosome sorting (numeric, then X, Y, MT)
        if chrom.replace('chr', '').isdigit():
            chrom_sort = int(chrom.replace('chr', ''))
        elif chrom.endswith('X'):
            chrom_sort = 23
        elif chrom.endswith('Y'):
            chrom_sort = 24
        elif chrom.endswith('M') or chrom.endswith('MT'):
            chrom_sort = 25
        else:
            chrom_sort = 26
        
        return (chrom_sort, pos)
    
    return sorted(merged_vars.items(), key=get_sort_key)


def main():
    parser = argparse.ArgumentParser(
        description='Merge VCF files from MuTect, LoFreqSomatic, and Strelka',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Priority rules:
  1. Variant in all 3 callers        -> Use MuTect INFO/FORMAT
  2. Variant in MuTect + Strelka     -> Use MuTect INFO/FORMAT
  3. Variant in LoFreq + Strelka     -> Use LoFreq INFO/FORMAT
  4. Variant in only 1 caller        -> Use that caller's INFO/FORMAT
        """
    )
    
    parser.add_argument('-m', '--mutect', required=True,
                        help='MuTect VCF file')
    parser.add_argument('-l', '--lofreqsomatic', required=True,
                        help='LoFreqSomatic VCF file')
    parser.add_argument('-s', '--strelka', required=True,
                        help='Strelka VCF file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output merged VCF file')
    parser.add_argument('--add-source-tag', action='store_true',
                        help='Add SOURCE tag to INFO field indicating caller(s)')
    
    args = parser.parse_args()
    
    # Read all VCF files
    print("Reading VCF files...", file=sys.stderr)
    mutect_header, mutect_vars = read_vcf(args.mutect, 'MuTect')
    lofreqsomatic_header, lofreqsomatic_vars = read_vcf(args.lofreqsomatic, 'LoFreqSomatic')
    strelka_header, strelka_vars = read_vcf(args.strelka, 'Strelka')
    
    # Merge headers
    print("Merging headers...", file=sys.stderr)
    merged_header = merge_headers(mutect_header, lofreqsomatic_header, strelka_header)
    
    # Get expected number of samples from merged header
    expected_samples = get_sample_count(merged_header)
    
    # Add SOURCE INFO field if requested
    if args.add_source_tag:
        source_line = '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Variant caller(s) that detected this variant">'
        # Insert before #CHROM line
        merged_header.insert(-1, source_line)
    
    # Merge variants
    print("Merging variants...", file=sys.stderr)
    merged_vars = merge_variants(mutect_vars, lofreqsomatic_vars, strelka_vars)
    
    # Sort variants
    print("Sorting variants...", file=sys.stderr)
    sorted_vars = sort_variants(merged_vars)
    
    # Write output
    print(f"Writing {len(sorted_vars)} variants to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as out:
        # Write header
        for line in merged_header:
            out.write(line + '\n')
        
        # Write variants
        for var_key, (fields, source) in sorted_vars:
            # Normalize variant line to have correct number of sample columns
            fields = normalize_variant_line(fields, expected_samples)
            
            if args.add_source_tag:
                # Add SOURCE to INFO field
                info_field = fields[7]
                if info_field == '.':
                    fields[7] = f'SOURCE={source}'
                else:
                    fields[7] = f'{info_field};SOURCE={source}'
            
            out.write('\t'.join(fields) + '\n')
    
    print(f"Done! Merged VCF written to {args.output}", file=sys.stderr)
    
    # Print summary statistics
    print("\nSummary:", file=sys.stderr)
    sources = defaultdict(int)
    for _, (_, source) in sorted_vars:
        sources[source] += 1
    
    for source in sorted(sources.keys()):
        print(f"  {source}: {sources[source]} variants", file=sys.stderr)


if __name__ == '__main__':
    main()