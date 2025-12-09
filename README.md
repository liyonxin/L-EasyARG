# L-EasyARG: A Comprehensive Pipeline for ARG Analysis from Metagenomic Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**L-EasyARG** is a robust and user-friendly bioinformatics pipeline designed for the identification and analysis of Antibiotic Resistance Genes (ARGs), Mobile Genetic Elements (MGEs), and Plasmids from metagenomic long-read sequencing data (e.g., Oxford Nanopore Technologies). It integrates taxonomic classification to provide a holistic view of the resistome.

## 🚀 Features

- **Comprehensive Profiling**: Simultaneous identification of ARGs (CARD), Plasmids (PLSDB), and MGEs.
- **Taxonomic Integration**: Links resistance genes to their host taxa using Centrifuge.
- **Automated Pipeline**: Streamlined workflow from raw reads to final merged results.
- **Visualization**: Built-in plotting capabilities for result interpretation.
- **Reproducibility**: Configuration-based execution for consistent results.

## 🛠️ Installation

### Prerequisites

Ensure the following external tools are installed and available in your system PATH:
- `seqkit`
- `seqtk`
- `centrifuge`
- `minimap2`
- `last` (lastal, lastdb, last-train)
- `maf-convert`
- `Rscript` (for some plotting functions)

### Install via Pip

You can install L-EasyARG directly from the source:

```bash
git clone https://github.com/liyonxin/L-EasyARG.git
cd L-EasyARG
pip install .
```

## 📚 Database Setup

L-EasyARG requires several reference databases. You can set them up automatically using the `setup` command:

```bash
easy-arg setup --database-dir /path/to/databases
```

This will download and configure:
- **CARD**: Comprehensive Antibiotic Resistance Database
- **PLSDB**: Plasmid Database
- **Centrifuge**: Bacterial/Viral/Archaeal index

> **Note**: The MGE database (mobileOG-db) may need to be downloaded manually if the automatic link is deprecated. Please place it in `L-EasyARG-database/mges/mobile-OG/`.

## 💻 Usage

### 1. Initialize Configuration

Create a default configuration file to customize parameters (threads, thresholds, database paths):

```bash
easy-arg init --output config.yaml
```

### 2. Run Analysis

Run the complete pipeline on your input data (FASTQ format):

```bash
easy-arg run --input sample.fastq.gz --output results_dir --threads 16
```

Or use a custom configuration file:

```bash
easy-arg run --input sample.fastq.gz --config config.yaml
```

### 3. Visualization

Generate plots from the analysis results:

```bash
easy-arg plot --input-dir results_dir --output-dir results_dir/plots --type all
```

## 📂 Output Structure

The pipeline generates the following directory structure:

- `rawdata/`: Symlinks/copies of input data and stats.
- `centrifuge/`: Taxonomic classification results.
- `ARG/`: ARG alignment results (PAF format).
- `plsdb/`: Plasmid alignment results.
- `MGE/`: Mobile Genetic Element alignments.
- `merged/`: **Final merged table** linking ARGs, Plasmids, MGEs, and Taxonomy.
- `plots/`: Generated visualizations.

## 📄 Citation

If you use L-EasyARG in your research, please cite:

> Li Y, et al. "L-EasyARG: A comprehensive tool for analyzing antibiotic resistance genes from metagenomic data." 

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
