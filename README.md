# Automated Workflow (Snakemake)

An optional `Snakefile` is provided to coordinate the variant consistency analysis steps and manage software dependencies through Snakemake.

This workflow serves as a convenience layer on top of the core analysis pipeline and may be used in place of manual execution of the individual scripts.

---

## 1. Environment Setup

The automated workflow requires a dedicated Snakemake controller environment to run the Snakemake engine.

**Hoffman2 Users:**  
Request an interactive compute node before creating environments to avoid login node timeouts.

Create and activate the Snakemake environment:

```
conda env create -f envs/snakemake.yaml
conda activate snakemake_env
```

---

## 2. Configuration

Before running the workflow, define your samples and global parameters in `config.yaml`.

---

### Sample Mapping

The workflow does not assume any specific naming convention for raw data files.

In the `samples:` section of `config.yaml`, provide:

- A unique ID for each sample
- The full path to its corresponding:
  - BAM file
  - barcodes file

**Example structure:**

```yaml
samples:
  sample_01:
    bam: /path/to/sample_01.bam
    barcodes: /path/to/sample_01/barcodes.tsv.gz

  sample_02:
    bam: /path/to/sample_02.bam
    barcodes: /path/to/sample_02/barcodes.tsv.gz
```

Each sample ID will be used to generate a corresponding output directory.

---

### Global Parameters

In addition to sample definitions, specify the following global parameters:

| Parameter | Description |
|------------|------------|
| `vcf` | Path to the reference genotype VCF |
| `donors` | Path to the `donors.txt` file containing multiplexed donor IDs |
| `out_root` | Directory where all numerical output folders (`00–03`) will be created |
| `coverage_thresholds` | List of coverage values used for final metric generation (e.g., `[0, 10, 20]`) |

Example:

```yaml
vcf: /path/to/reference.vcf.gz
donors: /path/to/donors.txt
out_root: results/
coverage_thresholds: [0, 10, 20]
```

---

## 3. Execution

Launch the pipeline with:

```
snakemake --use-conda --cores 16 -p
```

On the first run, Snakemake will automatically create the required sub-environments:

- `cellsnp-env`: For running `cellsnp-lite` during step 0
- `var_con-env`: For running the variant consistency analysis scripts during steps 1–3

---

## 4. Workflow Structure

The automated workflow executes the following staged directory structure:

```
00_cellsnp/   # Raw variant pileups from cellsnp-lite
01_counts/    # Generation of AD and DP matrices
02_indices/   # Partitioning of variants into C1, C2, I1, I2 categories
03_metrics/   # Final consistency CSVs (I1 variants used as proxy for ambient contamination)
```

This staged structure mirrors the manual execution steps while ensuring reproducible and scalable processing across samples and parameter settings.
