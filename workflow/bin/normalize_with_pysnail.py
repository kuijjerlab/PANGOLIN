import argparse
import os
from pysnail import Dataset, qsmooth

def parse_args():
    parser = argparse.ArgumentParser(description="Normalize expression data using qsmooth")
    parser.add_argument("xprs", type=str, help="Path to expression count matrix")
    parser.add_argument("groups", type=str, help="Path to sample groups file")
    parser.add_argument("output", type=str, help="Path to output normalized count file")
    parser.add_argument("--threshold", type=float, default=0.2, help="Threshold for qsmooth (default: 0.2)")
    return parser.parse_args()

def main():
    args = parse_args()

    xprs_path = os.path.realpath(args.xprs)
    groups_path = os.path.realpath(args.groups)
    output_path = os.path.realpath(args.output)

    dataset = Dataset(xprs_path, groups_path, index_col=0, sep='\t')
    xprs_norm, qstat = qsmooth(dataset, aggregation='auto', threshold=args.threshold)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    xprs_norm.to_csv(output_path, sep='\t')

if __name__ == "__main__":
    main()


