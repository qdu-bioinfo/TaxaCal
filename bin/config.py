import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='TaxaCal')

    ## Data settings
    parser.add_argument('--train_16S_g', type=str, default='data/training_16s_genus.csv')
    parser.add_argument('--train_WGS_g', type=str, default='data/training_wgs_genus.csv')
    parser.add_argument('--train_WGS_s', type=str, default='data/training_wgs_species.csv')
    parser.add_argument('--test', type=str, default='data/test_16s.csv')
    parser.add_argument('--o', type=str, default='data/calibrated.csv')

    return parser.parse_args()

