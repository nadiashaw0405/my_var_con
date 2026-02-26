Variant Consistency Pipeline
This pipeline evaluates donor-specific variant consistency in multiplexed single-cell data by comparing observed pileups against a reference VCF. It automates the process from raw BAM pileups to final consistency metrics using Snakemake.

1. Setup & Installation
The pipeline uses Snakemake to manage software environments automatically.

Clone this repository:

Bash
git clone <your-fork-url>
cd my_var_con
Create the Master Environment:
If you are on Hoffman2, load the Anaconda module first:

Bash
module load anaconda3
conda env create -f envs/snakemake.yaml
Activate the Environment:

Bash
conda activate snakemake_env
2. Configuration
Before running, you must update config.yaml with your specific project paths:

vcf: Path to your reference VCF.

donors: Path to a .txt file containing your multiplexed donor IDs.

samples: A list of your sample names.

out_root: The base directory where all results will be saved.

3. Running the Pipeline
Launch the pipeline with the following command. Snakemake will automatically build the required sub-environments for cellsnp-lite and the Python scripts.

Bash
snakemake --use-conda --cores 16 -p
Use -n before running to perform a "dry run" and verify the execution plan.

4. Output Structure
The pipeline generates an organized, numerical directory tree for each sample within your out_root:

00_cellsnp/: Variant pileup results.

01_counts/: AD and DP matrices.

02_indices/: Partitioned variant index dictionaries.

03_metrics/cov{X}/: Final consistency CSVs filtered by coverage threshold X.
