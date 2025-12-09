import os
import subprocess
import sys
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def check_dependencies():
    """Check if required tools are installed"""
    required_tools = [
        "seqkit", "seqtk", "centrifuge", "minimap2", 
        "lastal", "lastdb", "maf-convert", "Rscript", "wget", "tar", "unzip"
    ]
    
    missing = []
    for tool in required_tools:
        # Use shutil.which for a more pythonic check if available, 
        # but subprocess 'which' is fine for now.
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            missing.append(tool)
    
    if missing:
        logger.error(f"Missing required tools: {', '.join(missing)}")
        logger.error("Please install them before running L-EasyARG")
        sys.exit(1)

def run_command(cmd, shell=True, capture_output=False):
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
        logger.error(f"Error running command: {cmd}")
        logger.error(f"Error message: {e.stderr}")
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
        logger.error(f"Error parsing PAF file: {e}")
    
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
        logger.error(f"Error parsing centrifuge file: {e}")
    
    return taxid_dict

def setup_databases(database_dir, skip_download=False):
    """Setup required databases"""
    logger.info(f"Setting up databases in {database_dir}")
    os.makedirs(database_dir, exist_ok=True)
    
    # Database URLs and commands
    # Updated URLs based on setup.sh and known sources
    databases = {
        "card": {
            "url": "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2",
            "cmd": "tar -xjf broadstreet-v4.0.1.tar.bz2 && mv broadstreet-v4.0.1 card_database",
            "check_file": "card_database/nucleotide_fasta_protein_homolog_model.fasta"
        },
        "plsdb": {
            "url": "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta",
            "cmd": "mv download_fasta sequences.fasta",
            "check_file": "sequences.fasta"
        },
        "centrifuge": {
            "url": "https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed+h+v.tar.gz",
            "cmd": "tar -xzf p_compressed+h+v.tar.gz",
            "check_file": "p_compressed+h+v.1.cf"
        },
        "mge": {
             # Note: MGE database URL was not explicit in setup.sh, assuming it needs to be provided or downloaded separately.
             # For now, we'll keep a placeholder or use a known one if available.
             # The config points to: L-EasyARG-database/mges/mobile-OG/mobileOG-db_beatrix-1.6.All.faa
             # This looks like a specific structure.
             # I will add a warning if MGE is not found, as I don't have a direct URL for it in the snippets.
             "url": None, 
             "cmd": "echo 'Please download MGE database manually'",
             "check_file": "mges/mobile-OG/mobileOG-db_beatrix-1.6.All.faa"
        }
    }
    
    for db_name, db_info in databases.items():
        db_path = os.path.join(database_dir, db_name)
        if db_name == "mge":
             # MGE structure is deeper
             db_path = database_dir # Just use base dir for now as the config expects deep path
        
        os.makedirs(db_path, exist_ok=True)
        
        logger.info(f"\nSetting up {db_name} database...")
        
        # Check if already exists
        check_path = os.path.join(database_dir, db_name, db_info['check_file']) if db_name != "mge" else os.path.join(database_dir, db_info['check_file'])
        
        if os.path.exists(check_path):
             logger.info(f"  {db_name} database already exists. Skipping.")
             continue

        if not skip_download:
            if db_info['url']:
                # Download
                logger.info(f"  Downloading from {db_info['url']}")
                try:
                    download_cmd = f"cd {db_path} && wget -q {db_info['url']}"
                    subprocess.run(download_cmd, shell=True, check=True)
                    
                    # Extract/setup
                    logger.info(f"  Extracting/setting up...")
                    setup_cmd = f"cd {db_path} && {db_info['cmd']}"
                    subprocess.run(setup_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to setup {db_name}: {e}")
            else:
                logger.warning(f"  No download URL for {db_name}. Please install manually.")
        else:
            logger.info(f"  Skipping download (using existing files)")
    
    logger.info("\n✓ Database setup completed!")