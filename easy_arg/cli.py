#!/usr/bin/env python3
"""
L-EasyARG: A comprehensive tool for analyzing antibiotic resistance genes from metagenomic data
"""

import os
import sys
import argparse
import yaml
from pathlib import Path
from easy_arg.analysis import run_analysis_pipeline
from easy_arg.plotting import run_plotting_pipeline

def load_config(config_file=None):
    """Load configuration from file or use defaults"""
    default_config = {
        "database": {
            "card": "L-EasyARG-database/card_database/nucleotide_fasta_protein_homolog_model.fasta",
            "plsdb": "L-EasyARG-database/plsdb/sequences.fasta",
            "mge": "L-EasyARG-database/mges/mobile-OG/mobileOG-db_beatrix-1.6.All.faa",
            "centrifuge": "L-EasyARG-database/centrifuge/p+h+v/p_compressed+h+v",
            "who_species": "L-EasyARG-database/2024-WHO-species.txt"
        },
        "threads": 56,
        "min_identity": 0.75,
        "min_coverage": 0.7
    }
    
    if config_file and os.path.exists(config_file):
        with open(config_file, 'r') as f:
            user_config = yaml.safe_load(f)
            # Merge with defaults
            for key in user_config:
                if key in default_config and isinstance(default_config[key], dict):
                    default_config[key].update(user_config[key])
                else:
                    default_config[key] = user_config[key]
    
    return default_config

def setup_directories():
    """Create necessary directories"""
    dirs = ["rawdata", "dehost", "ARG", "MGE", "plsdb", "centrifuge", "merged", "R", "plots"]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(
        description="L-EasyARG: Analysis of Antibiotic Resistance Genes from Metagenomic Data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete analysis pipeline
  easy-arg run --input data.fastq.gz --sample SRR123456
  
  # Run with custom config
  easy-arg run --input data.fastq.gz --config my_config.yaml
  
  # Only run plotting on existing results
  easy-arg plot --type all
  
  # Run specific plot
  easy-arg plot --type top10 --sample SRR123456
        """
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Run command
    run_parser = subparsers.add_parser("run", help="Run complete analysis pipeline")
    run_parser.add_argument("--input", "-i", required=True, help="Input FASTQ file or directory")
    run_parser.add_argument("--sample", "-s", help="Sample name (default: basename of input)")
    run_parser.add_argument("--output", "-o", default=".", help="Output directory")
    run_parser.add_argument("--config", "-c", help="Configuration file")
    run_parser.add_argument("--threads", "-t", type=int, help="Number of threads")
    run_parser.add_argument("--skip-dehost", action="store_true", help="Skip host removal step")
    run_parser.add_argument("--skip-centrifuge", action="store_true", help="Skip taxonomic classification")
    
    # Plot command
    plot_parser = subparsers.add_parser("plot", help="Generate plots from analysis results")
    plot_parser.add_argument("--type", choices=["all", "top10", "distribution", "network", "cooccurrence"], 
                           default="all", help="Type of plot to generate")
    plot_parser.add_argument("--sample", "-s", help="Sample name for individual plots")
    plot_parser.add_argument("--input-dir", "-i", default=".", help="Input directory with analysis results")
    plot_parser.add_argument("--output-dir", "-o", default="plots", help="Output directory for plots")
    plot_parser.add_argument("--config", "-c", help="Configuration file")
    
    # Setup command
    setup_parser = subparsers.add_parser("setup", help="Setup databases")
    setup_parser.add_argument("--database-dir", "-d", default="L-EasyARG-database", 
                            help="Directory to store databases")
    setup_parser.add_argument("--skip-download", action="store_true", help="Skip download, only configure")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Load configuration
    config = load_config(args.config)
    
    if args.command == "run":
        # Update config with command line arguments
        if args.threads:
            config["threads"] = args.threads
        
        # Determine sample name
        if not args.sample:
            args.sample = Path(args.input).stem.replace(".fastq.gz", "").replace(".fq.gz", "")
        
        # Setup directories
        os.chdir(args.output)
        setup_directories()
        
        # Run analysis
        run_analysis_pipeline(
            input_path=args.input,
            sample_name=args.sample,
            config=config,
            skip_dehost=args.skip_dehost,
            skip_centrifuge=args.skip_centrifuge
        )
        
    elif args.command == "plot":
        # Run plotting
        run_plotting_pipeline(
            plot_type=args.type,
            sample_name=args.sample,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            config=config
        )
        
    elif args.command == "setup":
        # Setup databases
        from easy_arg.utils import setup_databases
        setup_databases(args.database_dir, skip_download=args.skip_download)
    
    print("\n✓ Analysis completed successfully!")

if __name__ == "__main__":
    main()