import os
import subprocess
import sys
from pathlib import Path

def check_dependencies():
    """Check if required tools are installed"""
    required_tools = [
        "seqkit", "seqtk", "centrifuge", "minimap2", 
        "lastal", "lastdb", "maf-convert", "Rscript"
    ]
    
    missing = []
    for tool in required_tools:
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            missing.append(tool)
    
    if missing:
        print(f"Error: Missing required tools: {', '.join(missing)}")
        print("Please install them before running L-EasyARG")
        sys.exit(1)

def run_command(cmd, shell=False, capture_output=False):
    """Run shell command with error handling"""
    try:
        if shell:
            result = subprocess.run(cmd, shell=True, capture_output=capture_output, 
                                  text=True, check=True)
        else:
            result = subprocess.run(cmd.split(), capture_output=capture_output,
                                  text=True, check=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error message: {e.stderr}")
        raise

def parse_paf(paf_file, gene_type='ARG'):
    """Parse PAF format file"""
    hits = []
    try:
        with open(paf_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 11:
                    hits.append({
                        'read_id': parts[0],
                        'type': gene_type,
                        'gene_name': parts[5],
                        'gene_length': int(parts[6]),
                        'start': int(parts[2]),
                        'end': int(parts[3]),
                        'strand': parts[4]
                    })
    except Exception as e:
        print(f"Error parsing PAF file: {e}")
    
    return hits

def parse_centrifuge(centrifuge_file):
    """Parse centrifuge results"""
    taxid_dict = {}
    try:
        with open(centrifuge_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    read_id = parts[0]
                    taxid = parts[2]
                    taxid_dict[read_id] = taxid
    except Exception as e:
        print(f"Error parsing centrifuge file: {e}")
    
    return taxid_dict

def setup_databases(database_dir, skip_download=False):
    """Setup required databases"""
    print(f"Setting up databases in {database_dir}")
    os.makedirs(database_dir, exist_ok=True)
    
    # Database URLs and commands
    databases = {
        "card": {
            "url": "https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2",
            "cmd": "tar -xjf broadstreet-v3.2.8.tar.bz2 && mv broadstreet-v3.2.8 card_database"
        },
        "plsdb": {
            "url": "https://ccb-microbe.cs.uni-saarland.de/plsdb/plsdb.fna.gz",
            "cmd": "gunzip -c plsdb.fna.gz > sequences.fasta"
        },
        "centrifuge": {
            "url": "https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed+h+v.tar.gz",
            "cmd": "tar -xzf p_compressed+h+v.tar.gz"
        }
    }
    
    for db_name, db_info in databases.items():
        db_path = os.path.join(database_dir, db_name)
        os.makedirs(db_path, exist_ok=True)
        
        print(f"\nSetting up {db_name} database...")
        
        if not skip_download:
            # Download
            print(f"  Downloading from {db_info['url']}")
            download_cmd = f"cd {db_path} && wget -q {db_info['url']}"
            subprocess.run(download_cmd, shell=True, check=True)
            
            # Extract/setup
            print(f"  Extracting/setting up...")
            setup_cmd = f"cd {db_path} && {db_info['cmd']}"
            subprocess.run(setup_cmd, shell=True, check=True)
        else:
            print(f"  Skipping download (using existing files)")
    
    print("\n✓ Database setup completed!")