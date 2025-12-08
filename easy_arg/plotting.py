import os
import subprocess
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def run_plotting_pipeline(plot_type, sample_name=None, input_dir=".", output_dir="plots", config=None):
    """Run plotting pipeline"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    if plot_type in ["all", "top10"]:
        if sample_name:
            plot_top10_args(sample_name, input_dir, output_dir, config)
        else:
            print("Sample name required for top10 plot")
    
    if plot_type in ["all", "distribution"]:
        plot_arg_distribution(input_dir, output_dir, config)
    
    if plot_type in ["all", "network"]:
        plot_network(input_dir, output_dir, config)
    
    if plot_type in ["all", "cooccurrence"]:
        plot_cooccurrence(input_dir, output_dir, config)

def plot_top10_args(sample_name, input_dir, output_dir, config):
    """Plot top 10 ARG subtypes (Python implementation)"""
    print(f"Generating top10 plot for {sample_name}")
    
    # Read merged results
    merged_file = os.path.join(input_dir, "merged", f"{sample_name}_merged_results.tsv")
    if not os.path.exists(merged_file):
        print(f"Warning: Merged file not found: {merged_file}")
        return
    
    # Read sequence length
    length_file = os.path.join(input_dir, "sum_length.txt")
    sample_bp = 0
    if os.path.exists(length_file):
        with open(length_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[0] == sample_name:
                    sample_bp = int(parts[1])
                    break
    
    if sample_bp == 0:
        print(f"Warning: No base pair length found for {sample_name}")
        return
    
    # Load data
    df = pd.read_csv(merged_file, sep='\t')
    arg_df = df[df['type'] == 'ARG']
    
    if len(arg_df) == 0:
        print(f"No ARG data found for {sample_name}")
        return
    
    # Calculate abundance
    abundance_factor = 1e9 / sample_bp
    arg_counts = arg_df['gene_name'].value_counts()
    arg_abundance = arg_counts * abundance_factor
    
    # Get top 10
    top10 = arg_abundance.head(10)
    
    # Plot
    plt.figure(figsize=(12, 8))
    colors = ["#8968CD", "#71C671", "#7EC0EE", "#6B8E23", "#436EEE", 
              "#388E8E", "#218868", "#27408B", "#7CCD7C", "#CD853F"]
    
    bars = plt.barh(range(len(top10)), top10.values, color=colors, alpha=0.9)
    
    # Add labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        plt.text(width + max(top10.values)*0.01, bar.get_y() + bar.get_height()/2,
                f'{width:.2f}', ha='left', va='center', fontweight='bold')
    
    plt.title(f'Top 10 ARG subtypes - {sample_name}', fontsize=14)
    plt.xlabel('Abundance (copies/Gb)', fontsize=12)
    plt.ylabel('ARG Subtype', fontsize=12)
    plt.yticks(range(len(top10)), top10.index)
    plt.grid(axis='x', alpha=0.3)
    
    # Save
    output_file = os.path.join(output_dir, f"{sample_name}_top10_args.pdf")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved top10 plot: {output_file}")

def plot_arg_distribution(input_dir, output_dir, config):
    """Plot ARG distribution across WHO priority pathogens"""
    print("Generating ARG distribution plot...")
    
    # Call R script
    r_script = Path(__file__).parent.parent / "scripts" / "R" / "R2_arg_distribution.R"
    if r_script.exists():
        cmd = [
            "Rscript", str(r_script),
            "--input", input_dir,
            "--output", output_dir,
            "--who-list", config['database']['who_species']
        ]
        subprocess.run(cmd, check=True)
    else:
        print("Warning: R script not found, using Python implementation")
        # Python implementation here...

def plot_network(input_dir, output_dir, config):
    """Plot ARG-microbe co-occurrence network"""
    print("Generating network plot...")
    
    r_script = Path(__file__).parent.parent / "scripts" / "R" / "R3_network.R"
    if r_script.exists():
        cmd = ["Rscript", str(r_script), "--input", input_dir, "--output", output_dir]
        subprocess.run(cmd, check=True)

def plot_cooccurrence(input_dir, output_dir, config):
    """Plot ARG-MGE co-occurrence"""
    print("Generating co-occurrence plot...")
    
    r_script = Path(__file__).parent.parent / "scripts" / "R" / "R4_cooccurrence.R"
    if r_script.exists():
        cmd = ["Rscript", str(r_script), "--input", input_dir, "--output", output_dir]
        subprocess.run(cmd, check=True)