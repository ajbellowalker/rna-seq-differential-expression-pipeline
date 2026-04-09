![FastQC](https://img.shields.io/badge/workflow-FastQC-blue)
![R](https://img.shields.io/badge/env-R-green)
![License](https://img.shields.io/badge/license-MIT-yellow)

## 🚀 Quick Start

```bash
# Clone repo
git clone https://github.com/ajbellowalker/rna-seq-differential-expression-pipeline
cd rna-seq-differential-expression-pipeline
```

# RNA-seq Differential Expression Pipeline

## 📌 Overview

This repository contains an end-to-end RNA-seq analysis pipeline for processing and analysing human transcriptomic data. The workflow integrates command-line bioinformatics tools with R-based statistical analysis to identify differentially expressed genes and perform functional enrichment analysis.

The pipeline is designed with a focus on **reproducibility, modularity, and data quality**, reflecting best practices in bioinformatics and genomic data engineering.

## ⚙️ Tools & Technologies

### Command-line tools

- **FastQC** – Quality control of raw sequencing reads  
- **Hisat2** – Alignment of reads to reference genome  
- **Samtools** – File conversion, sorting, and indexing  

### R / Bioconductor packages

- **DESeq2** – Differential expression analysis  
- **pheatmap** – Heatmap visualisation  
- **topGO** – Gene Ontology enrichment analysis  
- **org.Hs.eg.db** – Gene annotation database  

## 🔬 Workflow

The pipeline consists of the following steps:

1. **Quality Control**
   - Assess raw FASTQ files using FastQC  

2. **Read Alignment**
   - Align sequencing reads to the reference genome using Hisat2  

3. **File Processing**
   - Convert SAM to BAM format using Samtools  
   - Sort and index BAM files  

4. **Gene Expression Quantification**
   - Generate gene-level count data  

5. **Differential Expression Analysis**
   - Perform statistical analysis using DESeq2  
   - Identify significantly differentially expressed genes  

6. **Data Visualisation**
   - Principal Component Analysis (PCA)  
   - Volcano plots  
   - Heatmaps of significant genes  

7. **Functional Enrichment Analysis**
   - Perform Gene Ontology (GO) enrichment using topGO  
   - Identify biologically relevant pathways  

## 📊 Outputs

All outputs are stored in the `results/` directory:

### Tables (`results/tables/`)

- Differential expression results  
- Significant gene lists  
- Normalised count matrices  
- GO enrichment results  

### Visualisations (`results/plots/`)

- PCA plots  
- Volcano plots  
- Heatmaps  

## 📁 Repository Structure

├── data/ # Input data and sample metadata    
├── scripts/ # Pipeline scripts (Bash + R)   
├── results/ # Output files (tables and plots)   
├── docs/ # Workflow documentation   
└── assets/  

## ▶️ Pipeline Workflow

### 1. Quality Control

```bash
bash scripts/qc.sh
```

### 2. Alignment

```bash
bash scripts/alignment.sh
```

### 3. Processing

```bash
bash scripts/processing.sh
```

### 4. Differential Expression Analysis

```R
source("scripts/deseq2_analysis.R")
```

## 🧠 Key Features

- Modular and reproducible pipeline design
- Integration of command-line and R-based analysis
- Automated result generation and structured outputs
- Statistical filtering and multiple testing correction
- Functional interpretation through GO enrichment

## 🔮 Future Improvements

- Workflow automation using Nextflow
- Integration with cloud platforms (AWS)
- Extension to single-cell RNA-seq analysis
- Incorporation of workflow management and logging

## ⚠️ Notes

- Raw FASTQ files are not included 
- File paths adapted for portability

## 🧑‍💻 Author

Ayemere J. Bellowalker  
Bioinformatics | Microbiome | Computational Biology

GitHub: <https://github.com/ajbellowalker>

## 📄 License

MIT
