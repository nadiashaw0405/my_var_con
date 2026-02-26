configfile: "config.yaml"

# Global references from config
VCF = config["vcf"]
DONORS = config["donors"]
OUT = config["out_root"]

rule all:
    input:
        # Uses .keys() to iterate through sample IDs defined in config dictionary
        expand("{out}/{s}/03_metrics/cov{c}/{s}_c1_df_{c}.csv", 
               out=OUT, s=config["samples"].keys(), c=config["coverage_thresholds"])

rule cellsnp_00:
    input:
        # Pulls specific paths from config dictionary based on current sample
        bam = lambda wildcards: config["samples"][wildcards.sample]["bam"],
        barcodes = lambda wildcards: config["samples"][wildcards.sample]["barcodes"]
    output: directory("{out}/{sample}/00_cellsnp")
    log: "{out}/{sample}/logs/00_cellsnp.log"
    threads: 12
    conda: "envs/cellsnp.yaml"
    shell:
        "cellsnp-lite -s {input.bam} -b {input.barcodes} -R {VCF} -O {output} -p {threads} --gzip > {log} 2>&1"

rule run_script_01:
    input:  
        dir = "{out}/{sample}/00_cellsnp"
    output: 
        ad = "{out}/{sample}/01_counts/cellSNP.tag.AD.mtx.gz",
        dp = "{out}/{sample}/01_counts/cellSNP.tag.DP.mtx.gz"
    log: "{out}/{sample}/logs/01_counts.log"
    threads: 8
    conda: "envs/var_con.yaml"
    shell:
        """
        python 01_con_counts_multithread.py \
            -c {input.dir} -i {VCF} -d {DONORS} \
            -o {OUT}/{wildcards.sample}/01_counts -t {threads} > {log} 2>&1
        """

rule run_script_02:
    input:  
        # Explicitly requires the matrices from Rule 01 to ensure proper rule sequencing
        ad = "{out}/{sample}/01_counts/cellSNP.tag.AD.mtx.gz",
        dp = "{out}/{sample}/01_counts/cellSNP.tag.DP.mtx.gz"
    output: directory("{out}/{sample}/02_indices")
    log: "{out}/{sample}/logs/02_indices.log"
    threads: 16
    conda: "envs/var_con.yaml"
    shell:
        "python 02_get_con_indices.py -i {OUT}/{wildcards.sample}/01_counts -d {DONORS} -o {output} -t {threads} > {log} 2>&1"

rule run_script_03:
    input:  
        ad = "{out}/{sample}/01_counts/cellSNP.tag.AD.mtx.gz",
        idx_dir = "{out}/{sample}/02_indices"
    output: "{out}/{sample}/03_metrics/cov{cov}/{sample}_c1_df_{cov}.csv"
    log: "{out}/{sample}/logs/03_metrics_cov{cov}.log"
    conda: "envs/var_con.yaml"
    shell:
        """
        python 03_count_varcon_multithread.py \
            -i {OUT}/{wildcards.sample}/01_counts -p {input.idx_dir} -d {DONORS} \
            -c {wildcards.cov} \
            -o {OUT}/{wildcards.sample}/03_metrics/cov{wildcards.cov} > {log} 2>&1
        """
