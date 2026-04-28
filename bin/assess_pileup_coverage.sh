#!/usr/bin/env bash
# Median pileup depth (length of column 4 per row) vs minimum threshold.
# Mirrors the coverage gate in main.sh (median of mpileup-style depth per marker).
#
# Usage: assess_pileup_coverage.sh -T <tumor.pileup> -N <normal.pileup> [-m MIN]
#
# Exit status: 0 = PASS (both medians >= MIN), 1 = FAIL (either median < MIN or error).
# Stdout: PASS or FAIL (single line) for logging/parse.

set -euo pipefail

usage() {
	echo "Usage: $0 -T <tumor_pileup> -N <normal_pileup> [-m MIN_DEPTH]" >&2
	echo "  Default MIN_DEPTH is 50." >&2
}

TUMOR_PILEUP=""
NORMAL_PILEUP=""
COVERAGE_MIN=50

while getopts "T:N:m:h" opt; do # get options from command line
	case "$opt" in
	T) TUMOR_PILEUP=$OPTARG ;;
	N) NORMAL_PILEUP=$OPTARG ;;
	m) COVERAGE_MIN=$OPTARG ;;
	h) # help print usage and exit with 0 status
		usage
		exit 0
		;;
	*) usage; exit 1 ;; # exit with 1 status if invalid option
	esac # end case
done

# check if tumor and normal pileup files are provided
if [[ -z "$TUMOR_PILEUP" || -z "$NORMAL_PILEUP" ]]; then
	usage
	exit 1
fi

# check if pileup files exist
if [[ ! -f "$TUMOR_PILEUP" || ! -f "$NORMAL_PILEUP" ]]; then
	echo "ERROR: pileup file(s) not found" >&2
	exit 1
fi

# Median of mpileup depth per row = length of column 4 (read bases)
median_pileup_depth() {
	awk 'NF >= 4 { print length($4) }' "$1" | sort -n | awk '
		{ a[NR] = $1 }
		END {
			# if no rows, print 0 and exit
			if (NR == 0) { print 0; exit }
			# if odd number of rows, print the middle row
			if (NR % 2 == 1) print a[(NR + 1) / 2]
			# if even number of rows, print the average of the two middle rows
			else print (a[NR / 2] + a[NR / 2 + 1]) / 2
		}'
}

# calculate median pileup depth for tumor and normal samples
tumor_med=$(median_pileup_depth "$TUMOR_PILEUP")
normal_med=$(median_pileup_depth "$NORMAL_PILEUP")

# check if median pileup depth for tumor and normal samples is less than coverage minimum
if awk -v t="$tumor_med" -v n="$normal_med" -v m="$COVERAGE_MIN" 'BEGIN { exit !((t < m) || (n < m)) }'; then
	echo "FAIL"
	exit 1
fi

echo "PASS"
exit 0
