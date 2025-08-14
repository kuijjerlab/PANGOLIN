import argparse
import os
import pandas as pd
import random
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
import gc

def parse_args():
    parser = argparse.ArgumentParser(description="Run PANDA and LIONESS network inference")
    parser.add_argument("--exp_file", type=str, help="Path to expression matrix file")
    parser.add_argument("--motif_file", type=str, help="Path to motif file")
    parser.add_argument("--ppi_file", type=str, help="Path to PPI (protein-protein interaction) file")
    parser.add_argument("--output_dir", type=str, help="Output directory for results")
    parser.add_argument("--panda_output", type=str, default="panda_net.txt", 
                       help="PANDA output filename (default: panda_net.txt)")
    parser.add_argument("--start_sample", type=int, default=1, 
                       help="Start sample index for LIONESS (default: 1)")
    parser.add_argument("--end_sample", type=int, default=None, 
                       help="End sample index for LIONESS (default: all samples)")
    parser.add_argument("--computing", type=str, default="cpu", choices=["cpu", "gpu"],
                       help="Computing mode: cpu or gpu (default: cpu)")
    parser.add_argument("--ncores", type=int, default=1,
                       help="Number of cores to use (default: 1)"), # this is only for cpu
    parser.add_argument("--save_memory", action="store_true", 
                       help="Enable memory saving mode")
    parser.add_argument("--random_seed", type=int, default=10, 
                       help="Random seed for reproducibility (default: 10)")
    parser.add_argument("--lioness_precision", type=int, default=3, 
                       help="Precision for LIONESS output (default: 3)")
    parser.add_argument("--skip_lioness", action="store_true", 
                       help="Skip LIONESS analysis, only run PANDA")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Set random seed for reproducibility
    random.seed(args.random_seed)
    
    # Convert paths to absolute paths
    exp_file = os.path.realpath(args.exp_file)
    motif_file = os.path.realpath(args.motif_file)
    ppi_file = os.path.realpath(args.ppi_file)
    output_dir = os.path.realpath(args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read expression data
    print(f"Reading expression data from: {exp_file}")
    exp_df = pd.read_csv(exp_file, sep='\t', header=None, index_col=0)
    print(f"Expression matrix shape: {exp_df.shape}")
    
    # Run PANDA
    print("Running PANDA network inference...")
    panda_obj = Panda(exp_df, motif_file, ppi_file, 
                      save_tmp=False, 
                      save_memory=args.save_memory, 
                      remove_missing=False,
                      keep_expression_matrix=True, 
                      computing=args.computing)
    
    # Save PANDA results
    panda_output_path = os.path.join(output_dir, args.panda_output)
    print(f"Saving PANDA results to: {panda_output_path}")
    panda_obj.save_panda_results(panda_output_path)
    
    # Run LIONESS if not skipped
    if not args.skip_lioness:
        # Determine sample range
        if args.end_sample is None:
            num_samples = exp_df.shape[1]
            end_sample = num_samples
        else:
            end_sample = args.end_sample
        
        start_sample = args.start_sample
        
        print(f"Running LIONESS for samples {start_sample} to {end_sample}...")
        
        # Run LIONESS for each sample
        for sample_idx in range(start_sample, end_sample + 1):
            print(f"Processing sample {sample_idx}/{end_sample}")
            lioness_obj = Lioness(panda_obj, 
                                computing=args.computing, 
                                start=sample_idx, 
                                end=sample_idx, 
                                save_single=True, 
                                save_dir=output_dir,
                                save_fmt='txt', 
                                save_precision=args.lioness_precision, 
                                export_filename=None)
            
            # Clean up memory
            del lioness_obj
            gc.collect()
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main()
