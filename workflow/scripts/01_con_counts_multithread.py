#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack
from scipy.io import mmread, mmwrite

import subprocess

import argparse
import concurrent.futures
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Program to generate \
        sparse matrices containing consistency counts \
        per individual, for an input cell x SNP matrix (from cellSNP)')

parser.add_argument('-c', '--cellsnp_dir', type=str, required=True, \
        help='Pileup directory - direct output of cellSNP')

parser.add_argument('-i','--vcf', type=str, required=True, \
        help='VCF containing reference genotypes')

parser.add_argument('-d', '--donors', type=str, required=True, \
        help='txt file containing list of multiplexed donors')

parser.add_argument('-o', '--outdir', type=str, required=True, \
        help='directory to store sparse matrices per individual')

parser.add_argument('-t', '--threads', help='Number of threads', default=1)

args = parser.parse_args()

cellsnp_dir = args.cellsnp_dir
vcf_path = args.vcf
donor_path = args.donors
outdir = args.outdir
threads = int(args.threads)

###### function(s) ######

### breaks genotypes into hom/het
### based on agreement with pileup - consistency counts
def get_con_counts(ad, dp, geno):
    ### 0|1, 1|0
    if '0' in geno and '1' in geno:
        return dp
    ### 0|0
    elif '0' in geno and '1' not in geno:
        return dp - ad
    ### 1|1
    elif '0' not in geno and '1' in geno:
        return ad
    ### handles missing genotypes
    ### as well as multiallelic SNPs
    zeros = csr_matrix(np.zeros(dp.shape[1]))
    return zeros

def process_snp(ad, dp, genotypes):
  return get_con_counts(ad, dp, genotypes)

###### main I/O ######

donors = list(np.loadtxt(donor_path, dtype=str))

barcodes = pd.read_csv(f'{cellsnp_dir}/cellSNP.samples.tsv',
                       sep='\t',header=None,index_col=0)

# how do I detect this automatically?
demux_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
              'QUAL', 'FILTER', 'INFO']
# geno_base_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
#                   'QUAL', 'FILTER', 'INFO', 'FORMAT']
# geno_cols = geno_base_cols + donors

tmp_header = subprocess.run(f"zcat {vcf_path} | grep '#' | tail -1",
        shell=True, capture_output=True, text=True)
tmp_string = tmp_header.stdout
geno_cols = tmp_string.replace('\n', '').split('\t')

#
demux_snps = pd.read_csv(f'{cellsnp_dir}/cellSNP.base.vcf.gz',
      sep='\t', index_col=None, header=None, comment='#')
demux_snps.index = [f'{x}_{y}' for x,y in \
                    zip(demux_snps[0],demux_snps[1])]
demux_snps.columns = demux_cols

genotypes = pd.read_csv(vcf_path, sep='\t',
            index_col=None, header=None, comment='#')
genotypes.index = [f'{x}_{y}' for x,y in \
                   zip(genotypes[0],genotypes[1])]
genotypes.columns = geno_cols

# read sparse matrices
ad = mmread(f'{cellsnp_dir}/cellSNP.tag.AD.mtx').tocsr()
dp = mmread(f'{cellsnp_dir}/cellSNP.tag.DP.mtx').tocsr()
# TO DO: add functionality to process OTH alleles
oth = mmread(f'{cellsnp_dir}/cellSNP.tag.OTH.mtx').tocsr()

###### ensure that genotypes and cellSNP VCF use the same SNP set ######

d_mask = demux_snps.index.isin(genotypes.index)
g_mask = genotypes.index.isin(demux_snps.index)

tmp_genos = genotypes[g_mask]
tmp_demux = demux_snps[d_mask]

ad = ad[d_mask]
dp = dp[d_mask]
oth = oth[d_mask]

g_dup_mask = tmp_genos.index.duplicated()
d_dup_mask = tmp_demux.index.duplicated()

### probably easier way to get final ad/dp/oth but this works for now
ad = ad[~d_dup_mask]
dp = dp[~d_dup_mask]
oth = oth[~d_dup_mask]

final_genos = tmp_genos[~g_dup_mask]
final_demux = tmp_demux[~d_dup_mask]

assert final_genos.shape[0] == final_demux.shape[0]
assert dp.shape[0] == final_demux.shape[0]

n_snps = dp.shape[0]
ad_rows = [ad.getrow(i) for i in range(n_snps)]
dp_rows = [dp.getrow(i) for i in range(n_snps)]

for donor in donors:
  donor_genotypes = [x.split(':')[0] for x in final_genos[donor]]
  with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    mtx_list = list(tqdm(executor.map(process_snp,
                                      ad_rows, dp_rows, donor_genotypes),
                                      total=n_snps))
    consistent_mtx = vstack(mtx_list)
    mmwrite(f'{outdir}/{donor}.consistent.mtx', consistent_mtx)
    subprocess.run([f'gzip -f {outdir}/{donor}.consistent.mtx'], shell=True)

barcodes.to_csv(f'{outdir}/barcodes.tsv.gz', sep='\t',
        header=None, index=True)

mmwrite(f'{outdir}/cellSNP.tag.AD.mtx', ad)
mmwrite(f'{outdir}/cellSNP.tag.DP.mtx', dp)
mmwrite(f'{outdir}/cellSNP.tag.OTH.mtx', oth)

subprocess.run([f'gzip -f {outdir}/cellSNP.tag.AD.mtx'], shell=True)
subprocess.run([f'gzip -f {outdir}/cellSNP.tag.DP.mtx'], shell=True)
subprocess.run([f'gzip -f {outdir}/cellSNP.tag.OTH.mtx'], shell=True)

final_demux = final_demux.reset_index(drop=True)
final_demux.to_csv(f'{outdir}/varcon.SNPs.vcf.gz',
       sep='\t', header=True, index=True)
