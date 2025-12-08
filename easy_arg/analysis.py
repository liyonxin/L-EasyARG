import os
import subprocess
import pandas as pd
from pathlib import Path
from .utils import check_dependencies, run_command, parse_paf, parse_centrifuge

class AnalysisPipeline:
    def __init__(self, config):
        self.config = config
        self.path = os.getcwd()
        
    def run(self, input_path, sample_name, skip_dehost=False, skip_centrifuge=False):
        """Run complete analysis pipeline"""
        
        # Check dependencies
        check_dependencies()
        
        # Step 0: Prepare data
        self.prepare_data(input_path, sample_name)
        
        # Step 1: Run centrifuge (taxonomic classification)
        if not skip_centrifuge:
            self.run_centrifuge(sample_name)
        
        # Step 2: ARG identification
        self.run_arg_identification(sample_name)
        
        # Step 3: Plasmid identification
        self.run_plasmid_identification(sample_name)
        
        # Step 4: MGE identification
        self.run_mge_identification(sample_name)
        
        # Step 5: Filter results
        self.filter_results(sample_name)
        
        # Step 6: Merge results
        self.merge_results(sample_name)
        
        return True
    
    def prepare_data(self, input_path, sample_name):
        """Prepare input data and calculate statistics"""
        print(f"Preparing data for sample: {sample_name}")
        
        # Create symlink if input is a file
        if os.path.isfile(input_path):
            target = f"rawdata/{sample_name}.fastq.gz"
            if not os.path.exists(target):
                os.symlink(os.path.abspath(input_path), target)
        
        # Calculate sequence statistics
        cmd = f"seqkit stats rawdata/{sample_name}.fastq.gz -T"
        result = run_command(cmd, capture_output=True)
        
        # Parse and save stats
        with open(f"sum_length.txt", 'w') as f:
            f.write(f"{sample_name}\t{result.stdout.split()[3]}\n")
        
        # Convert to FASTA if needed
        fa_file = f"rawdata/{sample_name}.fa"
        if not os.path.exists(fa_file):
            cmd = f"seqtk seq -A rawdata/{sample_name}.fastq.gz > {fa_file}"
            run_command(cmd)
    
    def run_centrifuge(self, sample_name):
        """Run centrifuge for taxonomic classification"""
        print("Running centrifuge...")
        cmd = f"""centrifuge -f -x {self.config['database']['centrifuge']} \
                -U rawdata/{sample_name}.fa \
                --report-file centrifuge/{sample_name}_report.tsv \
                -S centrifuge/{sample_name}_result.tsv \
                -p {self.config['threads']}"""
        run_command(cmd)
    
    def run_arg_identification(self, sample_name):
        """Identify ARGs using minimap2"""
        print("Identifying ARGs...")
        cmd = f"""minimap2 -x map-ont --secondary=no \
                -t {self.config['threads']} \
                {self.config['database']['card']} \
                rawdata/{sample_name}.fa > ARG/{sample_name}_ARG.paf"""
        run_command(cmd)
    
    def run_plasmid_identification(self, sample_name):
        """Identify plasmids using minimap2"""
        print("Identifying plasmids...")
        cmd = f"""minimap2 -x map-ont --secondary=no \
                -t {self.config['threads']} \
                {self.config['database']['plsdb']} \
                rawdata/{sample_name}.fa > plsdb/{sample_name}_plsdb.paf"""
        run_command(cmd)
    
    def run_mge_identification(self, sample_name):
        """Identify MGEs using LAST"""
        print("Identifying MGEs...")
        
        # Create database if not exists
        mge_db = self.config['database']['mge']
        if not os.path.exists(mge_db + ".bck"):
            cmd = f"lastdb -P{self.config['threads']} -q -c trandb {mge_db}"
            run_command(cmd)
        
        # Run LAST
        cmds = [
            f"last-train -P{self.config['threads']} --codon trandb rawdata/{sample_name}.fa > MGE/{sample_name}.train",
            f"lastal -P{self.config['threads']} -p MGE/{sample_name}.train -m100 -D1e9 -K1 trandb rawdata/{sample_name}.fa > MGE/{sample_name}.maf",
            f"maf-convert psl MGE/{sample_name}.maf > MGE/{sample_name}_alignments.psl"
        ]
        for cmd in cmds:
            run_command(cmd)
    
    def filter_results(self, sample_name):
        """Filter alignment results"""
        print("Filtering results...")
        
        # Filter ARG results
        self.filter_paf(
            f"ARG/{sample_name}_ARG.paf",
            f"ARG/{sample_name}_ARG_filtered.txt",
            min_identity=0.75,
            min_coverage=0.9
        )
        
        # Filter plasmid results
        self.filter_paf(
            f"plsdb/{sample_name}_plsdb.paf",
            f"plsdb/{sample_name}_plsdb_filtered.txt",
            min_identity=0.7,
            min_coverage=0.7
        )
        
        # Filter MGE results
        self.filter_mge(sample_name)
    
    def filter_paf(self, input_file, output_file, min_identity, min_coverage):
        """Filter PAF format alignments"""
        if not os.path.exists(input_file):
            return
        
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) >= 11:
                    # Calculate identity and coverage
                    matches = int(parts[9])
                    aln_len = int(parts[10])
                    target_len = int(parts[6])
                    
                    identity = matches / aln_len if aln_len > 0 else 0
                    coverage = aln_len / target_len if target_len > 0 else 0
                    
                    if identity >= min_identity and coverage >= min_coverage:
                        f_out.write(line)
    
    def filter_mge(self, sample_name):
        """Filter MGE results"""
        input_file = f"MGE/{sample_name}_alignments.psl"
        output_file = f"MGE/{sample_name}_filtered_hits.txt"
        
        if not os.path.exists(input_file):
            return
        
        # Parse PSL format
        col_names = [
            'matches', 'mismatches', 'rep_matches', 'N_count', 'q_gap_count',
            'q_gap_bases', 't_gap_count', 't_gap_bases', 'strand', 'q_name',
            'q_len', 'q_start', 'q_end', 't_name', 't_len', 't_start',
            't_end', 'block_count', 'block_sizes', 'q_starts', 't_starts'
        ]
        
        try:
            df = pd.read_csv(input_file, sep='\t', header=None, names=col_names)
            df['protein_align_length'] = df['q_end'] - df['q_start']
            df['coverage'] = df['protein_align_length'] / df['q_len']
            df['identity'] = df['matches'] / df['protein_align_length']
            
            # Filter
            df_filtered = df[(df['coverage'] > 0.7) & (df['identity'] > 0.7)]
            
            # Remove overlapping hits
            df_filtered = self.remove_overlapping_hits(df_filtered)
            
            # Save
            output_cols = ['matches', 'strand', 'q_name', 'q_len', 'q_start',
                          'q_end', 't_name', 't_len', 't_start', 't_end', 'block_count']
            df_filtered[output_cols].to_csv(output_file, sep='\t', index=False)
            
        except Exception as e:
            print(f"Error filtering MGE results: {e}")
    
    def remove_overlapping_hits(self, df):
        """Remove overlapping MGE hits on same read"""
        if df.empty:
            return df
        
        df = df.sort_values(['t_name', 't_start'])
        df['to_remove'] = False
        
        for i in range(len(df) - 1):
            if df.iloc[i]['t_name'] == df.iloc[i + 1]['t_name']:
                # Calculate overlap
                overlap = min(df.iloc[i]['t_end'], df.iloc[i + 1]['t_end']) - \
                         max(df.iloc[i]['t_start'], df.iloc[i + 1]['t_start'])
                min_len = min(df.iloc[i]['t_end'] - df.iloc[i]['t_start'],
                             df.iloc[i + 1]['t_end'] - df.iloc[i + 1]['t_start'])
                
                if overlap / min_len > 0.8:
                    # Remove lower identity hit
                    if df.iloc[i]['identity'] > df.iloc[i + 1]['identity']:
                        df.at[df.index[i + 1], 'to_remove'] = True
                    else:
                        df.at[df.index[i], 'to_remove'] = True
        
        return df[~df['to_remove']].copy()
    
    def merge_results(self, sample_name):
        """Merge all results into a single file"""
        print("Merging results...")
        
        merged_data = []
        
        # Parse ARG results
        arg_file = f"ARG/{sample_name}_ARG_filtered.txt"
        if os.path.exists(arg_file):
            arg_hits = parse_paf(arg_file, 'ARG')
            merged_data.extend(arg_hits)
        
        # Parse MGE results
        mge_file = f"MGE/{sample_name}_filtered_hits.txt"
        if os.path.exists(mge_file):
            mge_hits = self.parse_mge_results(mge_file)
            merged_data.extend(mge_hits)
        
        # Parse centrifuge results
        centrifuge_file = f"centrifuge/{sample_name}_result.tsv"
        taxid_dict = {}
        if os.path.exists(centrifuge_file):
            taxid_dict = parse_centrifuge(centrifuge_file)
        
        # Parse plasmid results
        plasmid_reads = set()
        plasmid_file = f"plsdb/{sample_name}_plsdb_filtered.txt"
        if os.path.exists(plasmid_file):
            with open(plasmid_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if parts:
                        plasmid_reads.add(parts[0])
        
        # Write merged results
        output_file = f"merged/{sample_name}_merged_results.tsv"
        with open(output_file, 'w') as f:
            f.write("read_id\tplasmid_match\ttype\tgene_name\tgene_length\tstart\tend\tstrand\ttaxID\n")
            
            for hit in merged_data:
                read_id = hit['read_id']
                plasmid_match = "plasmid" if read_id in plasmid_reads else "genome"
                taxid = taxid_dict.get(read_id, "0")
                
                f.write(f"{read_id}\t{plasmid_match}\t{hit['type']}\t{hit['gene_name']}\t"
                       f"{hit['gene_length']}\t{hit['start']}\t{hit['end']}\t"
                       f"{hit['strand']}\t{taxid}\n")
    
    def parse_mge_results(self, mge_file):
        """Parse filtered MGE results"""
        hits = []
        try:
            df = pd.read_csv(mge_file, sep='\t')
            for _, row in df.iterrows():
                hits.append({
                    'read_id': row['t_name'],
                    'type': 'MGE',
                    'gene_name': row['q_name'],
                    'gene_length': row['q_len'],
                    'start': row['t_start'],
                    'end': row['t_end'],
                    'strand': row['strand'][1] if len(row['strand']) > 1 else '+'
                })
        except Exception as e:
            print(f"Error parsing MGE results: {e}")
        
        return hits

def run_analysis_pipeline(input_path, sample_name, config, **kwargs):
    """Main function to run analysis pipeline"""
    pipeline = AnalysisPipeline(config)
    return pipeline.run(input_path, sample_name, **kwargs)