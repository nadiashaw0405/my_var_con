# Variant Consistency Snakemake Workflow

This workflow coordinates the genotyping and variant consistency analysis for single-cell multiome data. It automates the transition from raw BAM files to donor-specific consistency metrics.

This workflow serves as a convenience layer on top of the core analysis pipeline and may be used in place of manual execution of the individual scripts.

---

## 1. Repository Structure

After cloning, your project directory should look like this:

```
my_var_con/
├── workflow/
│   ├── Snakefile             # Core automation logic and rules
│   ├── envs/                 # Conda environments
│   └── scripts/              # Python analysis scripts
├── config.yaml               # Sequencing run and reference file paths
├── submit_snakemake.sh       # Cluster submission script to run pipeline
├── logs/cluster/             # Cluster error/output logs
├── var_con_output/           # Results directory
└── README.md
```
---

## 2. Setup

### Environments
  1. **Base environment:** The workflow requires a base Snakemake environment. Create it once:
  ```
  conda env create -f workflow/envs/snakemake.yaml
  ```
  2. **Sub-environments:** The workflow requires two additional environments located in `workflow/envs/`.   It is recommended to manually pre-build the sub-environments from their respective `.yaml` files before running the pipeline for the first time. This ensures all dependencies are correctly indexed by your package manager:
  ```
  # 1. Install Mamba in your base environment if you don't have it
  conda install -n base -c conda-forge mamba
  
  # 2. Pre-build the pipeline environments
  mamba env create -f workflow/envs/cellsnp_env.yaml
  mamba env create -f workflow/envs/var_con_env.yaml
  ```
  
### Configuration
Edit `config.yaml` to define your runs and global parameters.
- `runs`: A mapping of run names to their respective BAM and barcode paths (see example below).
- `vcf`: Path to the reference genotype VCF (e.g., reference_000.vcf.gz).
- `donors`: A `.txt` file containing the donor IDs present in the multiplexed sample (see example below).
- `out_root`: The root directory for all outputs (default: `var_con_output`).
- `out_prefix`: The output directory for a specific run of the pipeline.
- `coverage`: Specify a minimum depth level, allowing comparison across quality thresholds (optional; default = `0`).

**Example structure for `runs`:**

```yaml
# Example config.yaml

runs: # Sequencing runs
  run_name_01:
    bam: "/path/to/run_name_01.bam"
    barcodes: "/path/to/run_name_01/barcodes.tsv.gz"

vcf: "/path/to/reference.vcf.gz"
donors: "/path/to/donors.txt"
out_root: "var_con_output"
out_prefix: "run_v1_test"
coverage: 0
```

**Example format for `donors.txt`:**

```text
Donor_A1
Donor_A2
Donor_A3
```
---

## 3. Execution

### Dry Run
Before launching the workflow, perform a dry run to verify that all paths in `config.yaml` are valid and the execution plan is correct.

```
conda activate snakemake_env
snakemake -np
```

This command prints the planned jobs without executing them.

### Local/Interactive Test

To verify that your BAM paths and Conda environments are working before launching a full run, execute the first step of the pipeline on an interactive compute node (via `qrsh`):
```
snakemake --use-conda --cores 12 cellsnp_00
```
If the pipeline says the directory is locked, run: `snakemake --unlock`.

### Launch Full Workflow on Cluster

For large-scale, automated processing of all sequencing runs defined in `config.yaml`, use the provided submission script to manage job scheduling:
```
qsub submit_snakemake.sh
```
**Monitor Progress:**
- **Main Pipeline Log:** `tail -f logs/snakemake_main.log`
- **Rule Logs:** `tail -f var_con_output/{prefix}/{run_name}/logs/00_cellsnp.log`
- **Cluster Logs:** `ls -ltr logs/cluster/`

---

## 4. Pipeline Stages & Results

| Stage | Rule | Tool | Output Description |
| :--- | :--- | :--- | :--- |
| **00** | `cellsnp_00` | `cellsnp-lite` | Genotypes single cells at VCF sites; outputs sparse AD/DP matrices. |
| **01** | `run_script_01` | `01_con_counts_multithread.py` | Generates donor-specific consistency matrices. |
| **02** | `run_script_02` | `02_get_con_indices.py` | Generates variant category indices (C1, C2, I1, I2) for donor partitioning. |
| **03** | `run_script_03` | `03_count_varcon_multithread.py` | Calculates final metrics and exports CSV results. |

After execution, results are written to the directory specified by `out_root` using the following structure:

```
{out_root}/                           # Root directory ("out_root" defined in config.yaml)
└── {out_prefix}/                     # Pipeline run prefix defined in "out_prefix," e.g., "run_v1"
    └── {run_name}/                   # Sequencing run name, e.g., "20201212-RUN-A1"
        │
        ├── 00_cellsnp/               # Step 0: Genotyping (cellsnp-lite)
        │   ├── cellSNP.base.vcf.gz   # VCF of sites evaluated
        │   ├── cellSNP.samples.tsv   # List of cell barcodes processed
        │   ├── cellSNP.tag.AD.mtx    # Allelic depth matrix 
        │   ├── cellSNP.tag.DP.mtx    # Total depth matrix 
        │   └── cellSNP.tag.OTH.mtx   # Other/error allele counts
        │
        ├── 01_counts/                       # Step 1: Donor Partitioning
        │   ├── varcon.SNPs.vcf.gz           # VCF filtered for sites present in reference genotypes
        │   ├── barcodes.tsv.gz              # Barcodes matching the matrix indices
        │   ├── cellSNP.tag.AD.mtx.gz        # Matrices subset to reference SNPs
        │   ├── cellSNP.tag.DP.mtx.gz  
        │   ├── cellSNP.tag.OTH.mtx.gz   
        │   ├── {Donor_A}.consistent.mtx.gz  # Matrices of consistent SNPs
        │   └── {Donor_B}.consistent.mtx.gz  # (one file per donor in your list)
        │
        ├── 02_indices/               # Step 2: Category Mapping
        │   ├── c1_dict.pkl           # Consistent variants (target donor)
        │   ├── c2_dict.pkl           # Consistent variants (other donors)
        │   ├── i1_dict.pkl           # Inconsistent/ambient proxy 1
        │   └── i2_dict.pkl           # Inconsistent/ambient proxy 2
        │
        ├── 03_metrics/               # Step 3: Final Results
        │   ├── {run_name}_c1_df.csv  # Final consistency metrics
        │   ├── {run_name}_c2_df.csv    
        │   ├── {run_name}_i1_df.csv    
        │   └── {run_name}_i2_df.csv 
        │
        └── logs/                     # Rule-level logs
            ├── 00_cellsnp.log
            ├── 01_counts.log
            ├── 02_indices.log
            └── 03_metrics.log
```
