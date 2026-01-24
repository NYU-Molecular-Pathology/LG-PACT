#!/usr/bin/env python3
"""
Combine SNV and CNV VCF files into a single VCF file.

Features:
- Deduplicates and merges headers intelligently
- Handles different sample structures (tumor-only vs tumor-normal)
- Adds missing NORMAL column to CNV variants
- Supports both .vcf and .vcf.gz files automatically
- Sorts variants by chromosome and position

Usage:
    combine_vcf.py -s snv.vcf -c cnv.vcf -o combined.vcf
    combine_vcf.py -s snv.vcf.gz -c cnv.vcf.gz -o combined.vcf.gz
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


def smart_open(filename, mode='rt'):
    """Open file, automatically handling .gz compression."""
    return gzip.open(filename, mode) if filename.endswith('.gz') else open(filename, mode)


def parse_vcf(filename):
    """Parse VCF file in single pass. Returns (headers, column_line, variants)."""
    headers, column_line, variants = [], None, []
    
    with smart_open(filename) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('##'):
                headers.append(line)
            elif line.startswith('#CHROM'):
                column_line = line
            elif line:
                variants.append(line)
    
    return headers, column_line, variants


def chr_sort_key(chrom):
    """Convert chromosome to sortable key (handles chr1-22, X, Y, M)."""
    chrom = chrom[3:] if chrom.startswith('chr') else chrom
    
    if chrom.isdigit():
        return (0, int(chrom))
    elif chrom == 'X':
        return (1, 0)
    elif chrom == 'Y':
        return (1, 1)
    elif chrom in ('M', 'MT'):
        return (1, 2)
    else:
        return (2, chrom)


def sort_variants(variants):
    """Sort variants by chromosome and position."""
    def sort_key(line):
        fields = line.split('\t')
        return (chr_sort_key(fields[0]), int(fields[1]))
    
    return sorted(variants, key=sort_key)


def add_normal_column(variants):
    """Add ./. to NORMAL column for variants with only TUMOR data."""
    result = []
    for variant in variants:
        if len(variant.split('\t')) == 10:
            variant += '\t./.'
        result.append(variant)
    return result


def merge_headers(snv_headers, cnv_headers):
    """Merge and deduplicate VCF headers, resolving conflicts intelligently."""
    
    merged = {}
    headers_by_type = {
        'FORMAT': {}, 'INFO': {}, 'FILTER': {},
        'ALT': {}, 'contig': {}, 'other': {}
    }
    
    def process_header(line, source='SNV'):
        """Process single header line."""
        if line.startswith('##fileformat'):
            merged['fileformat'] = line
        elif line.startswith('##reference'):
            merged['reference'] = line
        elif line.startswith('##fileDate'):
            merged['fileDate'] = line
        elif line.startswith('##source'):
            merged[f'source_{source.lower()}'] = line
        elif line.startswith('##FORMAT=<ID='):
            match = re.search(r'ID=([^,>]+)', line)
            if match:
                fmt_id = match.group(1)
                if fmt_id == 'GQ' and 'Type=Integer' not in line and fmt_id in headers_by_type['FORMAT']:
                    return
                if fmt_id not in headers_by_type['FORMAT'] or 'Type=Integer' in line:
                    headers_by_type['FORMAT'][fmt_id] = line
        elif line.startswith('##INFO=<ID='):
            match = re.search(r'ID=([^,>]+)', line)
            if match and match.group(1) not in headers_by_type['INFO']:
                headers_by_type['INFO'][match.group(1)] = line
        elif line.startswith('##FILTER=<ID='):
            match = re.search(r'ID=([^,>]+)', line)
            if match and match.group(1) not in headers_by_type['FILTER']:
                headers_by_type['FILTER'][match.group(1)] = line
        elif line.startswith('##ALT=<ID='):
            match = re.search(r'ID=([^,>]+)', line)
            if match and match.group(1) not in headers_by_type['ALT']:
                headers_by_type['ALT'][match.group(1)] = line
        elif line.startswith('##contig=<ID='):
            match = re.search(r'ID=([^,>]+)', line)
            if match:
                contig_id = match.group(1)
                if contig_id not in headers_by_type['contig'] or 'length=' in line:
                    headers_by_type['contig'][contig_id] = line
        else:
            headers_by_type['other'][line] = line
    
    for line in snv_headers:
        process_header(line, 'SNV')
    for line in cnv_headers:
        process_header(line, 'CNV')
    
    result = []
    for key in ['fileformat', 'reference', 'fileDate', 'source_snv', 'source_cnv']:
        if key in merged:
            result.append(merged[key])
    
    for header_type in ['ALT', 'FILTER', 'FORMAT', 'INFO', 'contig', 'other']:
        result.extend(sorted(headers_by_type[header_type].values()))
    
    return result


def validate_files(*files):
    """Validate that input files exist."""
    missing = [f for f in files if not Path(f).exists()]
    if missing:
        for f in missing:
            print(f"Error: File not found: {f}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Combine SNV and CNV VCF files into a single merged VCF file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -s snv.vcf -c cnv.vcf -o combined.vcf
  %(prog)s -s snv.vcf.gz -c cnv.vcf.gz -o combined.vcf.gz
  %(prog)s -s snv.vcf -c cnv.vcf -o combined.vcf --no-sort
"""
    )
    
    parser.add_argument('-s', '--snv', required=True, metavar='FILE', help='SNV VCF file')
    parser.add_argument('-c', '--cnv', required=True, metavar='FILE', help='CNV VCF file')
    parser.add_argument('-o', '--output', required=True, metavar='FILE', help='Output VCF file')
    parser.add_argument('--no-sort', action='store_true', help='Keep original order')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s 3.0')
    
    args = parser.parse_args()
    
    # Validate inputs
    validate_files(args.snv, args.cnv)
    
    # Parse files
    snv_headers, snv_columns, snv_variants = parse_vcf(args.snv)
    cnv_headers, cnv_columns, cnv_variants = parse_vcf(args.cnv)
    
    # Merge headers
    headers = merge_headers(snv_headers, cnv_headers)
    
    # Process variants
    variants = snv_variants + cnv_variants
    variants = add_normal_column(variants)
    
    if not args.no_sort:
        variants = sort_variants(variants)
    
    # Write output
    try:
        with smart_open(args.output, 'wt') as f:
            for header in headers:
                f.write(f"{header}\n")
            f.write(f"{snv_columns}\n")
            for variant in variants:
                f.write(f"{variant}\n")
        
        if args.verbose:
            print(f"Created: {args.output}")
            print(f"  Total variants: {len(variants)} (SNV: {len(snv_variants)}, CNV: {len(cnv_variants)})")
        
    except IOError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()