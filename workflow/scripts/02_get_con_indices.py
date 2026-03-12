#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack, hstack
from scipy.io import mmread, mmwrite
from sklearn.preprocessing import binarize

import argparse
import concurrent.futures
from tqdm import tqdm
import pickle

parser = argparse.ArgumentParser(description='Program to partition \
        variants into consistency categories, \
        stored in one dictionary per category.')

parser.add_argument('-i', '--indir', type=str, required=True, \
        help='Input directory')

parser.add_argument('-d', '--donors', type=str, required=True, \
        help='txt file containing list of multiplexed donors')

parser.add_argument('-o', '--outdir', type=str, required=True, \
        help='directory to store sparse matrices per individual')

parser.add_argument('-t', '--threads', help='Number of threads', default=8)

args = parser.parse_args()

indir = args.indir
donor_path = args.donors
outdir = args.outdir
threads = int(args.threads)

###### function ######

### from /u/project-cluo/igvf/pilot/multiome/ipynb/var_consistency/improve_varcon_code/01b_get_snp_indices.ipynb
def con_indices_per_cell(bin_consistent, bin_inconsistent, snps, index):

    n_donors = len(bin_consistent)

    one_cell_bin_con = [mtx[:,index] for mtx in bin_consistent]
    merged_bin_con = hstack(one_cell_bin_con).todense()

    one_cell_bin_incon = [mtx[:,index] for mtx in bin_inconsistent]
    merged_bin_incon = hstack(one_cell_bin_incon).todense()

    c1_indices = []
    c2_indices = []
    i1_indices = []
    i2_indices = []

    ### C1 variants
    per_snp_con = np.ravel(np.sum(merged_bin_con, axis=1))
    c1_mask = per_snp_con == 1
    tmp_c1_indices = merged_bin_con[c1_mask]
    tmp_c1_snps = snps.loc[c1_mask]

    tmp_c1_indices = np.ravel(np.argmax(tmp_c1_indices, axis=1))
    for i in range(n_donors):
        tmp_c1_mask = np.where(tmp_c1_indices == i)[0]
        tmp_snps = tmp_c1_snps.iloc[tmp_c1_mask]
        c1_indices.append(tmp_snps.index)

    ### C2
    c2_mask = per_snp_con > 1
    tmp_c2_indices = merged_bin_con[c2_mask]
    tmp_c2_snps = snps.loc[c2_mask]

    for i in range(n_donors):
        tmp_c2_mask = np.ravel(tmp_c2_indices[:,i]) > 0
        tmp_snps = tmp_c2_snps.iloc[tmp_c2_mask]
        c2_indices.append(tmp_snps.index)

    ### I1
    per_snp_incon = np.ravel(np.sum(merged_bin_incon, axis=1))
    i1_mask = (per_snp_incon < n_donors) & (per_snp_incon > 0)
    tmp_i1_indices = merged_bin_incon[i1_mask]
    tmp_i1_snps = snps.loc[i1_mask]

    for i in range(n_donors):
        tmp_i1_mask = np.ravel(tmp_i1_indices[:,i]) == 1
        tmp_snps = tmp_i1_snps.iloc[tmp_i1_mask]
        i1_indices.append(tmp_snps.index)

    ### I2
    ### handle sites that have both REF/ALT and only homozygous genos

    i2_mask = per_snp_incon == n_donors

    tmp_i2_indices = merged_bin_incon[i2_mask]
    tmp_i2_snps = snps.loc[i2_mask]

    ### bonus
    i2_incon = merged_bin_incon[i2_mask]
    i2_con = merged_bin_con[i2_mask]
    ### binarized - 2 indicates it's both consistent and inconsistent
    con_sums = (i2_incon + i2_con) < 2

    for i in range(n_donors):
        tmp_indices = np.ravel(con_sums[:,i])
        tmp_snps = tmp_i2_snps.iloc[tmp_indices]
        i2_indices.append(tmp_snps.index)

    return c1_indices, c2_indices, i1_indices, i2_indices

def process_barcode(i):
  return con_indices_per_cell(bin_consistent, bin_inconsistent, snps, i)

###### main I/O ######

donors = list(np.loadtxt(donor_path, dtype=str))
n_donors = len(donors)
barcodes = pd.read_csv(f'{indir}/barcodes.tsv.gz',
                       sep='\t',header=None,index_col=0)

# dp = mmread(f'{indir}/cellSNP.tag.DP.mtx').tocsr()
# consistent = [mmread(f'{indir}/{donor}.consistent.mtx').tocsr()
#               for donor in donors]

dp = mmread(f'{indir}/cellSNP.tag.DP.mtx.gz').tocsr()
consistent = [mmread(f'{indir}/{donor}.consistent.mtx.gz').tocsr()
              for donor in donors]

inconsistent = [dp - mtx for mtx in consistent]

bin_consistent = [binarize(x) for x in consistent]
bin_inconsistent = [binarize(x) for x in inconsistent]

snps = pd.read_csv(f'{indir}/varcon.SNPs.vcf.gz', sep='\t', header=0, index_col=0)
snps['chrom_pos'] = snps.index
snps.reset_index(drop=True, inplace=True)

n_barcodes = barcodes.shape[0]

### get indices first

c1_indices = []
c2_indices = []
i1_indices = []
i2_indices = []

with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    results = list(tqdm(executor.map(process_barcode, range(n_barcodes)), total=n_barcodes))

for c1_idx, c2_idx, i1_idx, i2_idx in results:
    c1_indices.append(c1_idx)
    c2_indices.append(c2_idx)
    i1_indices.append(i1_idx)
    i2_indices.append(i2_idx)

c1_dict = dict()
c2_dict = dict()
i1_dict = dict()
i2_dict = dict()

for i in range(n_barcodes):
    bc = barcodes.index[i]
    for j in range(n_donors):
        donor = donors[j]
        bc_donor = f'{bc}_{donor}'
        c1_dict[bc_donor] = c1_indices[i][j]
        c2_dict[bc_donor] = c2_indices[i][j]
        i1_dict[bc_donor] = i1_indices[i][j]
        i2_dict[bc_donor] = i2_indices[i][j]

with open(f'{outdir}/c1_dict.pkl', 'wb') as f:
    pickle.dump(c1_dict, f)
with open(f'{outdir}/c2_dict.pkl', 'wb') as f:
    pickle.dump(c2_dict, f)
with open(f'{outdir}/i1_dict.pkl', 'wb') as f:
    pickle.dump(i1_dict, f)
with open(f'{outdir}/i2_dict.pkl', 'wb') as f:
    pickle.dump(i2_dict, f)
