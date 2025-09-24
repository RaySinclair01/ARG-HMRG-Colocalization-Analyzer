# ARG-HMRG-Colocalization-Analyzer

A comprehensive bioinformatics pipeline to identify and analyze the co-localization of Antibiotic Resistance Genes (ARGs) and Heavy Metal Resistance Genes (HMRGs) in metagenomic data.

This workflow takes assembled contigs as input and produces statistical reports, abundance tables, and formatted data ready for gene cluster visualization. The pipeline is structured as a series of modular Bash and Python scripts, ensuring reproducibility and flexibility.

---

## 🏢 Institution

This project was developed at the:  
*Hunan Provincial University Key Laboratory for Environmental and Ecological Health*  
*Hunan Provincial University Key Laboratory for Environmental Behavior and Control Principle of New Pollutants*  
*College of Environment and Resources, Xiangtan University, Xiangtan 411105, China*

---


---

## ⚙️ Requirements

### Software
*   [Prodigal](https://github.com/hyattpd/Prodigal) for gene prediction.
*   [DIAMOND](https://github.com/bbuchfink/diamond) for protein sequence alignment.
*   A Unix-like environment (Linux, macOS, WSL on Windows) with Bash.
*   [Python 3](https://www.python.org/) with the following libraries:
    *   `pandas`
    *   `numpy` (usually installed with pandas)

### Databases
*   **CARD**: The Comprehensive Antibiotic Resistance Database.
*   **BacMet**: A database of antibacterial biocide and metal resistance genes.

---

## 🚀 Usage Instructions

It is recommended to run the scripts in sequential order as they depend on the output of the previous steps.

### Initial Setup

1.  Place your assembled contig files (e.g., `A1A_contigs.fasta`, `A1B_contigs.fasta`) into an `input_data/` directory.
2.  Prepare the CARD and BacMet databases as required by DIAMOND (`diamond makedb`).

### Script Workflow

**07. `07_create_prodigal_script.sh`**
*   **Purpose**: Predicts protein-coding genes from contigs using Prodigal.
*   **Action**: Generates and executes a script to run Prodigal in parallel for all samples defined in the `SAMPLES` array.
*   **Output**: Predicted proteins (`.faa`), genes (`.fna`), and coordinates (`.gff`) in the `analysis_results/` directory.

**08. `08_create_database_script.sh`** (Run once)
*   **Purpose**: Downloads and formats the CARD and BacMet databases.
*   **Action**: Generates and executes a script to download, combine (for BacMet), and build DIAMOND databases.
*   **Output**: DIAMOND database files (`.dmnd`) in a `~/databases/` directory.

**09. `09_create_arg_annotation_script.sh`**
*   **Purpose**: Annotates predicted proteins against the CARD database to find ARGs.
*   **Action**: Generates and executes a script to run DIAMOND for all samples.
*   **Output**: ARG hit tables (`_card_hits.m8`) in `analysis_results/`.

**10. `10_create_mge_annotation_script.sh`**
*   **Purpose**: Annotates predicted proteins against the BacMet database to find HMRGs.
*   **Action**: Generates and executes a script to run DIAMOND for all samples.
*   **Output**: HMRG hit tables (`_bacmet_hits.tsv`) in `analysis_results/`.

**11. `11_run_colocalization_analysis_final.py`**
*   **Purpose**: The core data integration step. Combines gene coordinates (GFF) with ARG/HMRG annotations.
*   **Action**: Identifies co-localized contigs and generates a detailed report for each sample, including correctly parsed gene names.
*   **Output**: Detailed co-localization reports (`_colocalization_full_details.tsv`) in a `colocalization_analysis_final_v3/` directory.

**12. `12_analyze_abundance_final.py`**
*   **Purpose**: Performs statistical analysis on the integrated data.
*   **Action**: Calculates the abundance of each unique "ARG-HMRG" pair based on the number of contigs they share.
*   **Output**: Gene pair abundance reports (`_abundance_report.tsv`) in an `abundance_analysis_final/` directory.

**13. `13_format_data_for_plotting.py`**
*   **Purpose**: Prepares the **complete** co-localization data for visualization tools.
*   **Action**: Converts the detailed report from script #11 into a clean, five-column format (`ID`, `source`, `start`, `end`, `strand`), while retaining annotation columns for advanced plotting.
*   **Output**: Formatted data (`_plot_data.tsv`) in a `cluster_plotting_data_final/` directory.

**14. `14_filter_by_contig_density.py`**
*   **Purpose**: Filters the dataset to focus on the most significant co-localization events.
*   **Action**: Calculates a "co-localization density score" for each contig and selects the Top N (e.g., Top 10) highest-scoring contigs. Includes a smart de-duplication step for genes with dual annotations.
*   **Output**: A smaller, focused dataset (`_top_10_contigs_plot_data.tsv`) in a `top_10_contigs_for_plotting/` directory.

**15. `15_format_for_external_plotter.py`**
*   **Purpose**: Formats the filtered, high-density data into a generic format for external visualization websites or software.
*   **Action**: Converts the output from script #14 into the strict five-column format required by many plotters, renaming contig IDs to generic group names (Contig_1, Contig_2, etc.).
*   **Output**: The final, ready-to-upload data (`_final_plot_data.tsv`) in a `final_data_for_visualization/` directory.

---

## 📄 Citing

If you use this pipeline in your research, please consider citing the tools it relies on, such as Prodigal, DIAMOND, CARD, and BacMet.

## 📧 Contact

For questions or issues with the pipeline, please open an issue on this GitHub repository.
