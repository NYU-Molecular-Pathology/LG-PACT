#!/usr/bin/env python3
"""
Combined QC Parser and MySQL Upload (Direct Connection - No SSH)
Parses paired tumor/normal QC CSV and uploads to MySQL database via direct connection
"""

import csv
import pymysql
import pymysql.cursors
import sys
import argparse
import os
import re
from decimal import Decimal, InvalidOperation

# Shared missing-value set for database conversion helpers.
MISSING_VALUES = {"", "NA", "N/A", "NULL", "NONE", "NO_DATA", "NAN"}

# =====================================================
# CONFIGURATION
# =====================================================
try:
    from db_config import DB_CONFIG

    try:
        from db_config import COMBINED_TABLE_NAME as TABLE_NAME
    except ImportError:
        from db_config import TABLE_NAME
    print("✓ Using configuration from config.py")
except ImportError:
    print("⚠ config.py not found, using default configuration")
    DB_CONFIG = {
        "host": "localhost",
        "user": "your_username",
        "password": "your_password",
        "database": "qc_database",
        "port": 3306,
    }
    TABLE_NAME = "sampleqc_metrics"


# =====================================================
# PARSING FUNCTIONS
# =====================================================


def parse_combined_qc_csv(input_file):
    """
    Parse paired tumor/normal QC CSV file
    Returns list of combined records (one per case)
    """
    parsed_data = []

    try:
        with open(input_file, "r") as infile:
            reader = csv.DictReader(infile)
            rows = list(reader)

            # Group rows by TM_number (case_id)
            cases = {}
            for row in rows:
                tm_number = row.get("TM_Number", "").strip()
                sample_type = row.get("Sample_Type", "").strip()

                if tm_number not in cases:
                    cases[tm_number] = {}

                cases[tm_number][sample_type] = row

            # Process each case
            for tm_number, samples in cases.items():
                tumor_row = samples.get("Tumor", {})
                normal_row = samples.get("Normal", {})

                if not tumor_row and not normal_row:
                    continue

                def get_value(row, key):
                    value = row.get(key, "").strip()
                    return "NA" if value == "" else value

                def clean_percentage(value):
                    if value == "NA" or value == "":
                        return "NA"
                    return value.replace("%", "").strip()

                def clean_snp_overlap(value):
                    if value == "NA" or value == "":
                        return "NA"
                    if "%" in value:
                        return value.split("%")[0].strip()
                    return value.strip()

                combined_record = {
                    "run_id": get_value(tumor_row or normal_row, "Run_ID"),
                    "case_id": tm_number if tm_number else "NA",
                    "tumor_id": get_value(tumor_row, "Case_ID"),
                    "normal_id": get_value(normal_row, "Case_ID"),
                    "status": get_value(tumor_row or normal_row, "Case_Level_QC"),
                    "tumor_status": get_value(tumor_row, "Sample_Overall_QC"),
                    "normal_status": get_value(normal_row, "Sample_Overall_QC"),
                    "tumor_mapped_reads": get_value(tumor_row, "Mapped Reads QC_Value"),
                    "tumor_dedup_reads": get_value(
                        tumor_row, "Deduplicated Reads QC_Value"
                    ),
                    "tumor_pct_target_ge_50x": clean_percentage(
                        get_value(tumor_row, "Targets with <50 Coverage_Value")
                    ),
                    "tumor_mean_coverage": get_value(tumor_row, "Average_Coverage"),
                    "tumor_median_coverage": get_value(tumor_row, "Median_Coverage"),
                    "tumor_snp_overlap": clean_snp_overlap(
                        get_value(tumor_row, "Overlapped HOMO SNPs_Value")
                    ),
                    "tumor_concordance": get_value(tumor_row, "Concordance_Percent"),
                    "normal_mapped_reads": get_value(
                        normal_row, "Mapped Reads QC_Value"
                    ),
                    "normal_dedup_reads": get_value(
                        normal_row, "Deduplicated Reads QC_Value"
                    ),
                    "normal_pct_target_ge_50x": clean_percentage(
                        get_value(normal_row, "Targets with <50 Coverage_Value")
                    ),
                    "normal_mean_coverage": get_value(normal_row, "Average_Coverage"),
                    "normal_median_coverage": get_value(normal_row, "Median_Coverage"),
                    "normal_snp_overlap": clean_snp_overlap(
                        get_value(normal_row, "Overlapped HOMO SNPs_Value")
                    ),
                    "normal_concordance": get_value(normal_row, "Concordance_Percent"),
                    "tumor_pct_target_bases_250x": clean_percentage(
                        get_value(tumor_row, "PCT_TARGET_BASES_250X")
                    ),
                    "tumor_exome_coverage": get_value(
                        tumor_row, "Exome_Percent_Coverage"
                    ),
                    "normal_pct_target_bases_250x": clean_percentage(
                        get_value(normal_row, "PCT_TARGET_BASES_250X")
                    ),
                    "normal_exome_coverage": get_value(
                        normal_row, "Exome_Percent_Coverage"
                    ),
                }

                parsed_data.append(combined_record)

        print(f"✓ Parsed {len(parsed_data)} case(s) from {input_file}")
        return parsed_data

    except FileNotFoundError:
        print(f"✗ Error: File '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except KeyError as e:
        print(f"✗ Error: Required column not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error parsing file: {e}", file=sys.stderr)
        sys.exit(1)


def save_parsed_csv(parsed_data, input_file):
    """Save parsed data to CSV file for reference"""
    try:
        base_name = os.path.basename(input_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = f"{name_without_ext}_parsed.csv"

        if parsed_data:
            fieldnames = [
                "run_id",
                "case_id",
                "tumor_id",
                "normal_id",
                "status",
                "tumor_status",
                "normal_status",
                "tumor_mapped_reads",
                "tumor_dedup_reads",
                "tumor_pct_target_ge_50x",
                "tumor_mean_coverage",
                "tumor_median_coverage",
                "tumor_snp_overlap",
                "tumor_concordance",
                "normal_mapped_reads",
                "normal_dedup_reads",
                "normal_pct_target_ge_50x",
                "normal_mean_coverage",
                "normal_median_coverage",
                "normal_snp_overlap",
                "normal_concordance",
                "tumor_pct_target_bases_250x",
                "tumor_exome_coverage",
                "normal_pct_target_bases_250x",
                "normal_exome_coverage",
            ]

            with open(output_file, "w", newline="") as outfile:
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(parsed_data)

            print(f"✓ Saved parsed CSV to: {output_file}")
            return output_file

    except Exception as e:
        print(f"⚠ Warning: Could not save parsed CSV: {e}")
        return None


# =====================================================
# DATABASE FUNCTIONS
# =====================================================
def create_connection(config):
    """Create a database connection using PyMySQL"""
    try:
        connection = pymysql.connect(
            host=config["host"],
            user=config["user"],
            password=config["password"],
            database=config["database"],
            port=config.get("port", 3306),
            cursorclass=pymysql.cursors.DictCursor,
        )
        print(f"✓ Connected to MySQL Server")
        print(f"✓ Using database: {config['database']}")
        return connection

    except pymysql.MySQLError as e:
        print(f"✗ Error connecting to MySQL: {e}", file=sys.stderr)
        sys.exit(1)


def verify_table_exists(connection, table_name):
    """Verify table exists - exit if it doesn't"""
    try:
        with connection.cursor() as cursor:
            # cursor.execute(f"SHOW TABLES LIKE '{table_name}'")
            cursor.execute("SHOW TABLES LIKE %s", (table_name,))
            result = cursor.fetchone()

        if not result:
            print(f"✗ Error: Table '{table_name}' does not exist", file=sys.stderr)
            print("This script only works with existing tables.")
            sys.exit(1)

        print(f"✓ Table '{table_name}' exists")

    except pymysql.MySQLError as e:
        print(f"✗ Error checking table: {e}", file=sys.stderr)
        sys.exit(1)


# Safely quote MySQL table names before using them in SQL strings.
def quote_identifier(identifier):
    """
    Safely quote a MySQL identifier such as a table name.

    SQL parameters can protect values, but not table names. Since table_name is
    interpolated into SQL strings, restrict it to ordinary identifier characters.
    """
    if not re.fullmatch(r"[A-Za-z0-9_]+", identifier):
        raise ValueError(f"Unsafe SQL identifier: {identifier}")

    return f"`{identifier}`"


# Build the logical duplicate keys from parsed upload records.
def get_upload_keys(parsed_data):
    """
    Return the logical upload keys for this file.

    A duplicate should be defined as the same case in the same run, not merely
    the same case_id appearing anywhere in the table.
    """
    keys = set()

    for row in parsed_data:
        run_id = str(row.get("run_id", "")).strip()
        case_id = str(row.get("case_id", "")).strip()

        if case_id in {"", "NA"}:
            continue

        if run_id in {"", "NA"}:
            continue

        keys.add((run_id, case_id))

    return sorted(keys)


def check_for_duplicates(connection, table_name, parsed_data):
    """
    Check whether records already exist for the same run_id + case_id pairs.

    This avoids false positives where the same case_id exists from another run.
    """
    try:
        upload_keys = get_upload_keys(parsed_data)

        if not upload_keys:
            return []

        safe_table_name = quote_identifier(table_name)
        placeholders = ", ".join(["(%s, %s)"] * len(upload_keys))
        params = []

        for run_id, case_id in upload_keys:
            params.extend([run_id, case_id])

        query = f"""
            SELECT run_id, case_id, COUNT(*) AS count
            FROM {safe_table_name}
            WHERE (run_id, case_id) IN ({placeholders})
            GROUP BY run_id, case_id
        """

        with connection.cursor() as cursor:
            cursor.execute(query, params)
            rows = cursor.fetchall()

        return [(row["run_id"], row["case_id"], row["count"]) for row in rows]

    except Exception as e:
        print(f"⚠ Could not check for duplicates: {e}")
        return []


# Delete existing records by run_id + case_id instead of case_id alone
def delete_existing_records(connection, table_name, record_keys):
    """
    Delete existing rows for specific run_id + case_id pairs.

    This prevents update mode from deleting the same case_id from unrelated runs.
    """
    try:
        if not record_keys:
            print("✓ No existing records to delete")
            return

        safe_table_name = quote_identifier(table_name)
        placeholders = ", ".join(["(%s, %s)"] * len(record_keys))
        params = []

        for run_id, case_id in record_keys:
            params.extend([run_id, case_id])

        query = f"""
            DELETE FROM {safe_table_name}
            WHERE (run_id, case_id) IN ({placeholders})
        """

        with connection.cursor() as cursor:
            cursor.execute(query, params)
            deleted = cursor.rowcount

        connection.commit()
        print(
            f"✓ Deleted {deleted} existing record(s) for {len(record_keys)} run/case pair(s)"
        )

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error deleting data: {e}", file=sys.stderr)
        sys.exit(1)


# Convert decimal-like percentage fields before database insert
def to_db_decimal(value):
    """
    Convert decimal-like values to Decimal, returning None for missing values.

    Handles values such as:
    98.7
    98.7%
    98.7 %
    1,234.5
    NA
    N/A
    NULL
    NO_DATA
    """
    if value is None:
        return None

    value = str(value).strip()

    if value.upper() in MISSING_VALUES:
        return None

    value = value.replace(",", "")
    value = value.replace("%", "")
    value = value.strip()

    try:
        return Decimal(value)

    except InvalidOperation as e:
        raise ValueError(f"Invalid decimal value for database insert: {value!r}") from e


# Convert integer-like fields safely before database insert.
def to_db_int(value):
    """Convert integer-like values to int, returning None for missing values."""
    if value is None:
        return None

    value = str(value).strip()

    if value.upper() in MISSING_VALUES:
        return None

    return int(float(value.replace(",", "")))


# EDITED: Convert float-like fields safely before database insert
def to_db_float(value):
    """Convert float-like values to float, returning None for missing values."""
    if value is None:
        return None

    value = str(value).strip()

    if value.upper() in MISSING_VALUES:
        return None

    value = value.replace(",", "")
    value = value.replace("%", "")
    value = value.strip()

    return float(value)


def insert_parsed_data(connection, parsed_data, table_name):
    """Insert parsed data into the database"""
    # EDITED: Quote table name before interpolating into the INSERT statement
    safe_table_name = quote_identifier(table_name)

    insert_query = f"""
    INSERT INTO {safe_table_name}
    (run_id, case_id, tumor_id, normal_id, status,
    tumor_status, normal_status,
    tumor_mapped_reads, tumor_dedup_reads, tumor_pct_target_ge_50x,
    tumor_mean_coverage, tumor_median_coverage, tumor_snp_overlap,
    tumor_concordance, tumor_pct_target_bases_250x, tumor_exome_coverage,
    normal_mapped_reads, normal_dedup_reads, normal_pct_target_ge_50x,
    normal_mean_coverage, normal_median_coverage, normal_snp_overlap,
    normal_concordance, normal_pct_target_bases_250x, normal_exome_coverage)
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    # EDITED: Convert missing text values to None; preserve valid text values
    def to_db_value(value):
        """Convert missing string values to None; preserve valid text values."""
        if value is None:
            return None

        value = str(value).strip()
        if value.upper() in MISSING_VALUES:
            return None
        return value

    try:
        records_inserted = 0

        with connection.cursor() as cursor:
            for row in parsed_data:
        # Concordance fields are converted to numeric decimals before insert
                values = (
                    to_db_value(row["run_id"]),
                    to_db_value(row["case_id"]),
                    to_db_value(row["tumor_id"]),
                    to_db_value(row["normal_id"]),
                    to_db_value(row["status"]),
                    to_db_value(row["tumor_status"]),
                    to_db_value(row["normal_status"]),
                    to_db_int(row["tumor_mapped_reads"]),
                    to_db_int(row["tumor_dedup_reads"]),
                    to_db_float(row["tumor_pct_target_ge_50x"]),
                    to_db_float(row["tumor_mean_coverage"]),
                    to_db_float(row["tumor_median_coverage"]),
                    to_db_float(row["tumor_snp_overlap"]),
                    to_db_decimal(row["tumor_concordance"]),
                    to_db_float(row["tumor_pct_target_bases_250x"]),
                    to_db_value(row["tumor_exome_coverage"]),
                    to_db_int(row["normal_mapped_reads"]),
                    to_db_int(row["normal_dedup_reads"]),
                    to_db_float(row["normal_pct_target_ge_50x"]),
                    to_db_float(row["normal_mean_coverage"]),
                    to_db_float(row["normal_median_coverage"]),
                    to_db_float(row["normal_snp_overlap"]),
                    to_db_decimal(row["normal_concordance"]),
                    to_db_float(row["normal_pct_target_bases_250x"]),
                    to_db_value(row["normal_exome_coverage"]),
                )

        # DEBUG: Prints the exact concordance values before database insert
                # print(
                #     "DEBUG concordance:",
                #     "run_id=",
                #     row["run_id"],
                #     "case_id=",
                #     row["case_id"],
                #     "tumor_concordance=",
                #     repr(row["tumor_concordance"]),
                #     "normal_concordance=",
                #     repr(row["normal_concordance"]),
                # )

                cursor.execute(insert_query, values)
                records_inserted += 1

        connection.commit()
        print(f"✓ Successfully inserted {records_inserted} records")

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error inserting data: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        connection.rollback()
        print(f"✗ Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


# Table stats no longer assumes the table has a created_at column
def get_table_stats(connection, table_name):
    """Get statistics about the table."""
    try:
        safe_table_name = quote_identifier(table_name)

        with connection.cursor() as cursor:
            cursor.execute(f"SELECT COUNT(*) AS total FROM {safe_table_name}")
            total = cursor.fetchone()["total"]

            cursor.execute(f"""
                SELECT status, COUNT(*) AS count
                FROM {safe_table_name}
                GROUP BY status
            """)
            status_counts = cursor.fetchall()

            cursor.execute(f"""
                SELECT run_id, COUNT(*) AS count
                FROM {safe_table_name}
                GROUP BY run_id
                ORDER BY run_id DESC
                LIMIT 5
            """)
            run_counts = cursor.fetchall()

        return {
            "total": total,
            "by_status": status_counts,
            "by_run": run_counts,
        }

    except pymysql.MySQLError as e:
        print(f"⚠ Could not get table stats: {e}")
        return None


# =====================================================
# MAIN FUNCTION
# =====================================================
def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description="Upload paired tumor/normal QC data to existing MySQL table",
        epilog="""
Modes:
  append (default)  - Add all records to existing data
  update            - Replace data for case_ids in CSV, keep other cases

Examples:
  python %(prog)s combined_qc.csv
  python %(prog)s combined_qc.csv --update
  python %(prog)s combined_qc.csv --table other_table
        """,
    )

    parser.add_argument("input_csv", help="Path to input combined QC CSV file")
    parser.add_argument(
        "--update",
        "-u",
        action="store_true",
        help="Update mode: delete existing data for case_ids in CSV before inserting",
    )
    parser.add_argument("--table", "-t", help="Table name (overrides config)")
    parser.add_argument(
        "--yes", "-y", action="store_true", help="Skip confirmation prompts"
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_csv):
        print(f"✗ Error: Input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)

    table_name = args.table if args.table else TABLE_NAME
    mode = "update" if args.update else "append"

    print("=" * 80)
    print("Combined QC Upload (Tumor/Normal Paired Data)")
    print("=" * 80)
    print(f"Input File: {args.input_csv}")
    print(f"Database:   {DB_CONFIG['database']}")
    print(f"Table:      {table_name}")
    print(f"Mode:       {mode}")
    print("=" * 80)

    # Step 1: Parse CSV --------------------------------------------------
    print("\n[Step 1/3] Parsing combined QC CSV file...")
    parsed_data = parse_combined_qc_csv(args.input_csv)

    if not parsed_data:
        print("✗ No data to upload")
        sys.exit(1)

    save_parsed_csv(parsed_data, args.input_csv)

    # Step 2: Connect and verify -------------------------------------------
    print("\n[Step 2/3] Connecting to database...")
    connection = create_connection(DB_CONFIG)

    try:
        verify_table_exists(connection, table_name)

    # updated to use the new duplicate check: considers run_id + case_id pairs
        print("\nChecking for existing data...")
        duplicates = check_for_duplicates(connection, table_name, parsed_data)

        if duplicates:
            print(f"⚠ Found existing data for {len(duplicates)} run/case pair(s):")
            for run_id, case_id, count in duplicates[:5]:
                print(f"  - run_id={run_id}, case_id={case_id}: {count} record(s)")

            if len(duplicates) > 5:
                print(f"  ... and {len(duplicates) - 5} more")
        else:
            print("✓ No duplicate run/case pairs found")

        # Step 3: Upload --------------------------------------------------
        print("\n[Step 3/3] Uploading data...")

    # Update/delete logic targets specific run_id + case_id pairs instead of just case_id
        if args.update and duplicates:
            record_keys_to_delete = [
                (run_id, case_id) for run_id, case_id, count in duplicates
            ]

            if not args.yes:
                print(
                    f"\n⚠ Update mode will delete existing data for "
                    f"{len(record_keys_to_delete)} run/case pair(s)"
                )
                response = input("Continue? Type 'yes' to confirm: ")

                if response.lower() != "yes":
                    print("Upload cancelled.")
                    sys.exit(0)

            delete_existing_records(connection, table_name, record_keys_to_delete)

        elif args.update and not duplicates:
            print("No existing run/case pairs to update, will append data")

        insert_parsed_data(connection, parsed_data, table_name)

        print("\n" + "=" * 80)
        print("TABLE STATISTICS")
        print("=" * 80)
        stats = get_table_stats(connection, table_name)
        if stats:
            print(f"\nTotal Records: {stats['total']}")

            if stats["by_status"]:
                print("\nRecords by Status:")
                for item in stats["by_status"]:
                    print(f"  {item['status']}: {item['count']}")

            if stats["by_run"]:
                print("\nRecent Run IDs (top 5):")
                for item in stats["by_run"]:
                    print(f"  {item['run_id']}: {item['count']} cases")

        print("\n" + "=" * 80)
        print("✓ Upload completed successfully!")
        print("=" * 80)

    finally:
        connection.close()
        print("\n✓ MySQL connection closed")


if __name__ == "__main__":
    main()
