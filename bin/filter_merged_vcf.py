#!/usr/bin/env python3
"""
Filter VCF for PASS variants with VAF >= threshold.
Handles AF in INFO field (LoFreq) and FORMAT field (MuTect).
"""

import sys
import argparse
import gzip


def open_vcf(filename):
    """Open VCF file (handles gzipped files)."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')


def collect_used_fields(variants):
    """
    Scan variants to find which INFO, FORMAT, and FILTER fields are actually used.
    
    Args:
        variants: List of variant lines (not including headers)
    
    Returns:
        Tuple of (used_info_fields, used_format_fields, used_filter_values)
    """
    used_info = set()
    used_format = set()
    used_filter = set()
    
    for line in variants:
        fields = line.rstrip().split('\t')
        
        # Collect FILTER values
        filter_val = fields[6]
        if filter_val != '.':
            # Handle multiple filters separated by semicolon
            for filt in filter_val.split(';'):
                used_filter.add(filt)
        
        # Collect INFO fields
        info_field = fields[7]
        if info_field != '.':
            for item in info_field.split(';'):
                if '=' in item:
                    field_name = item.split('=')[0]
                    used_info.add(field_name)
                else:
                    # Flag fields (no value)
                    used_info.add(item)
        
        # Collect FORMAT fields
        if len(fields) > 8:
            format_field = fields[8]
            if format_field != '.':
                for fmt in format_field.split(':'):
                    used_format.add(fmt)
    
    return used_info, used_format, used_filter


def filter_headers(header_lines, used_info, used_format, used_filter, filtering_applied):
    """
    Filter header lines to only include definitions for fields that are actually used.
    
    Args:
        header_lines: List of header lines from input VCF
        used_info: Set of INFO field names actually used
        used_format: Set of FORMAT field names actually used
        used_filter: Set of FILTER values actually used
        filtering_applied: String describing the filtering applied
    
    Returns:
        List of filtered header lines
    """
    filtered_headers = []
    
    for line in header_lines:
        # Always keep these
        if line.startswith('##fileformat'):
            filtered_headers.append(line)
            continue
        elif line.startswith('##reference'):
            filtered_headers.append(line)
            continue
        elif line.startswith('##contig'):
            filtered_headers.append(line)
            continue
        elif line.startswith('##ALT'):
            filtered_headers.append(line)
            continue
        elif line.startswith('#CHROM'):
            # Add filtering info before the column header line
            if filtering_applied:
                filtered_headers.append(f'##filter_vcf={filtering_applied}')
            filtered_headers.append(line)
            continue
        
        # Filter INFO lines
        elif line.startswith('##INFO=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_info:
                filtered_headers.append(line)
        
        # Filter FORMAT lines
        elif line.startswith('##FORMAT=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_format:
                filtered_headers.append(line)
        
        # Filter FILTER lines
        elif line.startswith('##FILTER=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_filter or field_id == 'PASS':  # Always keep PASS definition
                filtered_headers.append(line)
        
        # Keep other header lines
        else:
            filtered_headers.append(line)
    
    return filtered_headers


def is_strelka_only(info_field):
    """Check if variant is only called by Strelka based on SOURCE tag."""
    for item in info_field.split(';'):
        if item.startswith('SOURCE='):
            source = item.split('=')[1]
            return source == 'Strelka'
    return False  # If no SOURCE tag, don't filter out


def get_af_from_info(info_field):
    """Extract AF value from INFO field."""
    for item in info_field.split(';'):
        if item.startswith('AF='):
            return float(item.split('=')[1])
    return None


def get_af_from_format(format_field, sample_field):
    """Extract AF value from FORMAT field."""
    format_tags = format_field.split(':')
    sample_values = sample_field.split(':')
    
    if 'AF' in format_tags:
        af_index = format_tags.index('AF')
        if af_index < len(sample_values):
            try:
                return float(sample_values[af_index])
            except (ValueError, IndexError):
                return None
    return None


def get_variant_af(fields):
    """
    Get allele frequency from variant line.
    Checks INFO field first, then FORMAT fields.
    Returns the maximum AF found across all samples.
    """
    info_field = fields[7]
    
    # Try INFO field first (LoFreq style)
    af = get_af_from_info(info_field)
    if af is not None:
        return af
    
    # Try FORMAT fields (MuTect/Strelka style)
    if len(fields) > 9:
        format_field = fields[8]
        max_af = None
        
        # Check all sample columns (usually TUMOR is the one with AF > 0)
        for i in range(9, len(fields)):
            sample_af = get_af_from_format(format_field, fields[i])
            if sample_af is not None:
                if max_af is None or sample_af > max_af:
                    max_af = sample_af
        
        if max_af is not None:
            return max_af
    
    return None


def filter_vcf(input_file, output_file, min_af=0.01, filter_status='PASS', exclude_strelka_only=False, clean_headers=False):
    """
    Filter VCF file for variants passing criteria.
    
    Args:
        input_file: Input VCF filename
        output_file: Output VCF filename
        min_af: Minimum allele frequency (default 0.01 = 1%)
        filter_status: Filter status to keep (default 'PASS'). Use 'ALL' to keep all FILTER values.
        exclude_strelka_only: If True, exclude variants only called by Strelka (default False)
        clean_headers: If True, remove unused header definitions (default False)
    """
    kept = 0
    filtered_out = 0
    strelka_only_filtered = 0
    
    # Check if we should skip filter checking
    skip_filter = filter_status.upper() in ['ALL', '*', 'ANY']
    
    # Read all data first (needed for header cleaning)
    header_lines = []
    kept_variants = []
    
    with open_vcf(input_file) as infile:
        for line in infile:
            # Collect header lines
            if line.startswith('#'):
                header_lines.append(line.rstrip())
                continue
            
            # Parse variant line
            fields = line.rstrip().split('\t')
            
            # Check if Strelka-only and should be excluded
            if exclude_strelka_only:
                info_field = fields[7]
                if is_strelka_only(info_field):
                    strelka_only_filtered += 1
                    filtered_out += 1
                    continue
            
            # Check FILTER column (skip if user wants all)
            if not skip_filter:
                filter_col = fields[6]
                if filter_col != filter_status:
                    filtered_out += 1
                    continue
            
            # Check AF
            af = get_variant_af(fields)
            
            if af is None:
                # No AF found - keep the variant with a warning
                print(f"Warning: No AF found for {fields[0]}:{fields[1]}, keeping variant", 
                      file=sys.stderr)
                kept_variants.append('\t'.join(fields))
                kept += 1
            elif af >= min_af:
                kept_variants.append('\t'.join(fields))
                kept += 1
            else:
                filtered_out += 1
    
    # Clean headers if requested
    if clean_headers and kept_variants:
        used_info, used_format, used_filter = collect_used_fields(kept_variants)
        
        # Build filtering description
        filter_desc = f'min_AF={min_af}'
        if not skip_filter:
            filter_desc += f';FILTER={filter_status}'
        if exclude_strelka_only:
            filter_desc += ';exclude_Strelka_only=true'
        
        header_lines = filter_headers(header_lines, used_info, used_format, used_filter, filter_desc)
    
    # Write output
    with open(output_file, 'w') as outfile:
        # Write headers
        for line in header_lines:
            outfile.write(line + '\n')
        
        # Write variants
        for variant in kept_variants:
            outfile.write(variant + '\n')
    
    return kept, filtered_out, strelka_only_filtered


def main():
    parser = argparse.ArgumentParser(
        description='Filter VCF for PASS variants with VAF >= threshold',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete filtering: PASS + VAF >= 1% + no Strelka-only + clean headers
  python filter_vcf.py -i merged.vcf -o filtered.vcf -a 0.01 --exclude-strelka-only --clean-headers
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input VCF file (can be gzipped)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output VCF file')
    parser.add_argument('-a', '--min-af', type=float, default=0.01,
                        help='Minimum allele frequency (default: 0.01 = 1%%)')
    parser.add_argument('-f', '--filter', default='PASS',
                        help='Filter status to keep (default: PASS). Use "ALL" to keep all variants regardless of FILTER.')
    parser.add_argument('--exclude-strelka-only', action='store_true',
                        help='Exclude variants only called by Strelka (requires SOURCE tag in INFO field)')
    parser.add_argument('--clean-headers', action='store_true',
                        help='Remove unused INFO/FORMAT/FILTER definitions from headers')
    
    args = parser.parse_args()
    
    print(f"Filtering VCF...", file=sys.stderr)
    print(f"  Input: {args.input}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)
    print(f"  Min AF: {args.min_af} ({args.min_af*100}%)", file=sys.stderr)
    print(f"  Filter: {args.filter}", file=sys.stderr)
    if args.exclude_strelka_only:
        print(f"  Exclude Strelka-only: Yes", file=sys.stderr)
    if args.clean_headers:
        print(f"  Clean headers: Yes", file=sys.stderr)
    
    kept, filtered, strelka_filtered = filter_vcf(args.input, args.output, args.min_af, args.filter, 
                                                   args.exclude_strelka_only, args.clean_headers)
    
    print(f"\nResults:", file=sys.stderr)
    print(f"  Kept: {kept} variants", file=sys.stderr)
    print(f"  Filtered out: {filtered} variants", file=sys.stderr)
    if args.exclude_strelka_only and strelka_filtered > 0:
        print(f"    - Strelka-only: {strelka_filtered} variants", file=sys.stderr)
    print(f"  Total: {kept + filtered} variants", file=sys.stderr)
    print(f"\nDone! Filtered VCF written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
