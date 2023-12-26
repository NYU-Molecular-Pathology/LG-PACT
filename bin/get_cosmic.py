#!/usr/bin/env python3

import pandas as pd
import re

class AnnotationsProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.data = pd.read_csv(self.file_path, sep="\t")

    def extract_required_columns(self):
        self.data['Variant'] = self.data[['CHROM', 'POS', 'REF', 'ALT']].apply(lambda x: ':'.join(map(str, x)), axis=1)
        self.data = self.data[['CHROM', 'POS', 'REF', 'ALT', 'cosmic70', 'Variant']]
        self.data = self.data.astype({'POS': str})

    def filter_cosmic_entries(self):
        self.data = self.data[self.data['cosmic70'] != "."]

    def parse_cosmic_info(self):
        self.data[['CosmicID', 'Occurence']] = self.data['cosmic70'].str.split(';', expand=True)
        self.data['CosmicID'] = self.data['CosmicID'].str.split(',').str[0]
        self.data['CosmicID'] = self.data['CosmicID'].str.replace("ID=", "")
        self.data['all_occurences'] = self.data['Occurence'].str.findall(r'\d+')
        self.data['int_all_occurences'] = [list(map(int, val)) for val in self.data['all_occurences']]
        self.data['COSMIC_Count'] = [sum(int_vals) for int_vals in self.data['int_all_occurences']]
        self.data = self.data[["Variant","CosmicID","COSMIC_Count"]]

    def get_processed_data(self):
        return self.data

def main(annotations_per_caller_file):
    processor = AnnotationsProcessor(annotations_per_caller_file)
    processor.extract_required_columns()
    processor.filter_cosmic_entries()
    processor.parse_cosmic_info()
    processed_data = processor.get_processed_data()
    return processed_data 

if __name__ == "__main__":
    main()

