#!/usr/bin/env python

__author__ = "Kelsey Zhu"
__version__ = "1.0.1"
import os
import jinja2
import argparse
from create_seq_report import seq_qc_report

def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="NGS607 output directory")
    parser.add_argument("-r", "--run", type=str, required=True,
                        help="NGS RUN ID, e.g. 210715_NB501073_0191_AHHTHKBGXJ")
    parser.add_argument("-p", "--project", type=str, required=True,
                        help="NGS607 project name, e.g. PACT-21-15")
    parser.add_argument("-s", "--SeraCare", type=str, required=False,
                        help="snp overlaping percentage",
                        default="/gpfs/data/molecpathlab/ref/QC/SeraCare-selected-variants-dist.tsv")
    return parser.parse_args()

def main():
    args = get_options()
    qc_report = seq_qc_report(args).build_report()
    clinical_path = os.path.join(args.output,"output/clinical") 
    if not os.path.exists(clinical_path): os.mkdir(clinical_path,0o0755)

    template_path = os.path.join(args.output,"template/")
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "seq_report_template.html"

    template = templateEnv.get_template(TEMPLATE_FILE)
    outputText = template.render(run_df=qc_report['run_df'], df=qc_report['df'],
                                 seracare=qc_report['SeraCare'],run=qc_report['run'],caseID=qc_report["caseID"])
    html_file = open('%s/output/clinical/%s.html'%(args.output,qc_report["caseID"]), 'w')
    html_file.write(outputText)
    html_file.close()

if __name__ == "__main__":
    main()
