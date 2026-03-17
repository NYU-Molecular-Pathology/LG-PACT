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
        if line.startswith('##fileformat'):
            filtered_headers.append(line)
        elif line.startswith('##reference'):
            filtered_headers.append(line)
        elif line.startswith('##contig'):
            filtered_headers.append(line)
        elif line.startswith('##ALT'):
            filtered_headers.append(line)
        elif line.startswith('#CHROM'):
            if filtering_applied:
                filtered_headers.append(f'##filter_vcf={filtering_applied}')
            filtered_headers.append(line)
        elif line.startswith('##INFO=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_info:
                filtered_headers.append(line)
        elif line.startswith('##FORMAT=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_format:
                filtered_headers.append(line)
        elif line.startswith('##FILTER=<ID='):
            field_id = line.split('ID=')[1].split(',')[0]
            if field_id in used_filter or field_id == 'PASS':
                filtered_headers.append(line)
        else:
            filtered_headers.append(line)

    return filtered_headers


def get_tumor_col_index(header_lines):
    """
    Parse the #CHROM header line to find the index of the TUMOR sample column.

    The expected sample name is 'TUMOR' (case-insensitive match).
    Prints a warning and returns None if the column is not found.

    Args:
        header_lines: List of header lines from the VCF.

    Returns:
        Integer index into fields[] for the TUMOR sample column, or None.
    """
    for line in header_lines:
        if line.startswith('#CHROM'):
            cols = line.split('\t')
            # Sample columns start at index 9
            for i, name in enumerate(cols[9:], start=9):
                if name.strip().upper() == 'TUMOR':
                    return i
            print("Warning: No TUMOR sample column found in #CHROM header line. "
                  "AF from FORMAT fields will not be extracted.", file=sys.stderr)
            return None
    return None


def is_strelka_only(info_field):
    """Check if variant is only called by Strelka based on SOURCE tag."""
    for item in info_field.split(';'):
        if item.startswith('SOURCE='):
            source = item.split('=')[1]
            return source == 'Strelka'
    return False


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


def get_variant_af(fields, tumor_col):
    """
    Get allele frequency from variant line using only the TUMOR sample column.

    Checks INFO field first (LoFreq style — single-sample, no normal column),
    then falls back to the TUMOR sample FORMAT column (MuTect/Strelka style).
    The NORMAL sample column is intentionally ignored.

    Args:
        fields: Split fields of a VCF data line.
        tumor_col: Index into fields[] for the TUMOR sample column.

    Returns:
        Float AF value, or None if not found.
    """
    # Try INFO field first (LoFreq style)
    af = get_af_from_info(fields[7])
    if af is not None:
        return af

    # Try TUMOR sample FORMAT column only
    if tumor_col is not None and len(fields) > tumor_col:
        return get_af_from_format(fields[8], fields[tumor_col])

    return None


def filter_vcf(input_file, output_file, min_af=0.01, filter_status='PASS',
               exclude_strelka_only=False, clean_headers=False):
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

    skip_filter = filter_status.upper() in ['ALL', '*', 'ANY']

    header_lines = []
    kept_variants = []
    tumor_col = None  # resolved once after all headers are read

    with open_vcf(input_file) as infile:
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line.rstrip())
                continue

            # Resolve TUMOR column index once, after all header lines are collected
            if tumor_col is None:
                tumor_col = get_tumor_col_index(header_lines)
                print(f"  TUMOR sample column index: {tumor_col}", file=sys.stderr)

            fields = line.rstrip().split('\t')

            # Exclude Strelka-only variants if requested
            if exclude_strelka_only and is_strelka_only(fields[7]):
                strelka_only_filtered += 1
                filtered_out += 1
                continue

            # Check FILTER column
            if not skip_filter and fields[6] != filter_status:
                filtered_out += 1
                continue

            # Check AF (TUMOR sample only)
            af = get_variant_af(fields, tumor_col)

            if af is None:
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

        filter_desc = f'min_AF={min_af}'
        if not skip_filter:
            filter_desc += f';FILTER={filter_status}'
        if exclude_strelka_only:
            filter_desc += ';exclude_Strelka_only=true'

        header_lines = filter_headers(header_lines, used_info, used_format, used_filter, filter_desc)

    # Write output
    with open(output_file, 'w') as outfile:
        for line in header_lines:
            outfile.write(line + '\n')
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
                        help='Minimum allele frequency (default: 0.01 = 1%)')
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

    kept, filtered, strelka_filtered = filter_vcf(
        args.input, args.output, args.min_af, args.filter,
        args.exclude_strelka_only, args.clean_headers
    )

    print(f"\nResults:", file=sys.stderr)
    print(f"  Kept: {kept} variants", file=sys.stderr)
    print(f"  Filtered out: {filtered} variants", file=sys.stderr)
    if args.exclude_strelka_only and strelka_filtered > 0:
        print(f"    - Strelka-only: {strelka_filtered} variants", file=sys.stderr)
    print(f"  Total: {kept + filtered} variants", file=sys.stderr)
    print(f"\nDone! Filtered VCF written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()