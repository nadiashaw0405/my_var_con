configfile: "config.yaml"

# Global references
VCF = config["vcf"]
DONORS = config["donors"]
OUT = config["out_root"]

rule all:
    input:
        expand("{out}/{s}/03_metrics/cov{c}/{s}_c1_df_{c}.csv", 
               out=OUT, s=config["samples"], c=config["coverage_thresholds"])

rule cellsnp_00:
    input:
        bam = config["bam_path"],
        barcodes = config["barcodes_path"]
    output: directory("{out}/{sample}/00_cellsnp")
    log: "{out}/{sample}/logs/00_cellsnp.log"
    threads: 12
    conda: "envs/cellsnp.yaml"  # <--- Activates cellsnp-env
    shell:
        "cellsnp-lite -s {input.bam} -b {input.barcodes} -R {VCF} -O {output} -p {threads} --gzip > {log} 2>&1"

rule run_script_01:
    input:  dir = "{out}/{sample}/00_cellsnp"
    output: ad = "{out}/{sample}/01_counts/cellSNP.tag.AD.mtx.gz",
            dp = "{out}/{sample}/01_counts/cellSNP.tag.DP.mtx.gz"
    log: "{out}/{sample}/logs/01_counts.log"
    threads: 8
    conda: "envs/var_con.yaml"   # <--- Activates var_con-env
    shell:
        """
        python 01_con_counts_multithread.py \
            -c {input.dir} -i {VCF} -d {DONORS} \
            -o {OUT}/{wildcards.sample}/01_counts -t {threads} > {log} 2>&1
        """

rule run_script_02:
    input:  dir = "{out}/{sample}/01_counts"
    output: directory("{out}/{sample}/02_indices")
    log: "{out}/{sample}/logs/02_indices.log"
    threads: 16
    conda: "envs/var_con.yaml"   # <--- Activates var_con-env
    shell:
        "python 02_get_con_indices.py -i {input.dir} -d {DONORS} -o {output} -t {threads} > {log} 2>&1"

rule run_script_03:
    input:  c_dir = "{out}/{sample}/01_counts",
            idx_dir = "{out}/{sample}/02_indices"
    output: "{out}/{sample}/03_metrics/cov{cov}/{sample}_c1_df_{cov}.csv"
    log: "{out}/{sample}/logs/03_metrics_cov{cov}.log"
    conda: "envs/var_con.yaml"   # <--- Activates var_con-env
    shell:
        """
        python 03_count_varcon_multithread.py \
            -i {input.c_dir} -p {input.idx_dir} -d {DONORS} \
            -c {wildcards.cov} \
            -o {OUT}/{wildcards.sample}/03_metrics/cov{wildcards.cov} > {log} 2>&1
        """
