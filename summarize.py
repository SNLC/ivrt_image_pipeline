#!/usr/bin/env python

import os
import pandas as pd
import glob
import argparse


def parse_args():

    parser = argparse.ArgumentParser(
            description="""Combine individual coverslip quantificiation
                        data into an experiment summary""",
                        argument_default=None)

    parser.add_argument('--path',
                        type=str,
                        default=os.getcwd(),
                        help='path to folder containing CSVs')

    args = parser.parse_args()

    return args


def join_df(csv_list):
    df_list = []
    for csv in csv_list:
        df = pd.read_csv(csv)
        df.set_index(keys=['Starter Cell Type',
                           'Target Cell Type',
                           'Experiment ID',
                           'Coverslip ID'],
                     inplace=True)
        df_list.append(df)
    concatenated_summary_df = pd.concat(df_list)
    return concatenated_summary_df


def main():
    args = parse_args()
    args.path = '/Volumes/My Book/rabies_tracing_images/pv_cre_starter_cells/XH_12_07_17/vglut_presynaptic_647'
    summary_list = glob.glob(os.path.join(args.path, '**/summary.csv'),
                             recursive=True)
    summary_df = join_df(summary_list)

    summary_outpath = os.path.join(args.path, 'summarized_counts/summary.csv')
    summary_df.to_csv(summary_outpath)

if __name__ == '__main__':
    main()
