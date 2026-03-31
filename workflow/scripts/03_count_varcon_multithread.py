import numpy as np
import pandas as pd
import gzip
import pickle
import argparse
from scipy.io import mmread
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Program to take \
        dictionaries of variants/consistency counts and a coverage threshold \
        and returns a dataframe with all types of consistency (to add: \
        filter by other criteria, such as imputation R^2/coverage).')

parser.add_argument('-i', '--indir', type=str, required=True, \
        help='Input directory')

parser.add_argument('-p', '--pkldir', type=str, required=True, \
        help='Directory with consistency dicts')

parser.add_argument('-d', '--donors', type=str, required=True, \
        help='txt file containing list of multiplexed donors')

parser.add_argument('-c', '--cov_thresh', type=str, required=True, \
        help='experiment-wide coverage threshold')

parser.add_argument('-o', '--outdir', type=str, required=True, \
        help='directory to store sparse matrices per individual')

parser.add_argument('-n', '--run_name', type=str, required=True, 
        help='prefix for output files (e.g., the sequencing run ID)')

args = parser.parse_args()

indir = args.indir
pkldir = args.pkldir
donor_path = args.donors
cov_thresh = int(args.cov_thresh)
outdir = args.outdir
run_name = args.run_name

### multithreading function

def process_bc(i):
    tmp_c1 = c1_indices[i]
    tmp_c2 = c2_indices[i]
    tmp_i1 = i1_indices[i]
    tmp_i2 = i2_indices[i]

    tmp_c1_counts = []
    tmp_c2_counts = []
    tmp_i1_counts = []
    tmp_i2_counts = []

    for j in range(n_donors):
        c1 = vcf.iloc[tmp_c1[j]]
        c2 = vcf.iloc[tmp_c2[j]]
        i1 = vcf.iloc[tmp_i1[j]]
        i2 = vcf.iloc[tmp_i2[j]]

        c1_mask = c1['DP'] > cov_thresh
        c2_mask = c2['DP'] > cov_thresh
        i1_mask = i1['DP'] > cov_thresh
        i2_mask = i2['DP'] > cov_thresh

        tmp_c1_counts.append(np.sum(consistent[j][c1.index[c1_mask], i]))
        tmp_c2_counts.append(np.sum(consistent[j][c2.index[c2_mask], i]))

        i1_consistent_sum = np.sum(consistent[j][i1.index[i1_mask], i])
        i1_total_dp_sum = np.sum(dp[i1.index[i1_mask]])
        tmp_i1_counts.append(i1_total_dp_sum - i1_consistent_sum)

        i2_consistent_sum = np.sum(consistent[j][i2.index[i2_mask], i])
        i2_total_dp_sum = np.sum(dp[i2.index[i2_mask]])
        tmp_i2_counts.append(i2_total_dp_sum - i2_consistent_sum)

    return tmp_c1_counts, tmp_c2_counts, tmp_i1_counts, tmp_i2_counts


###### filter SNPs ######
### to do: incorporate other filters, such as R^2

def filter_snps(snps):
  pass

###### main I/O ######

donors = list(np.loadtxt(donor_path, dtype=str))
n_donors = len(donors)
barcodes = pd.Index(pd.read_csv(f'{indir}/barcodes.tsv.gz',
                       sep='\t',header=None,index_col=0).index)

header_line = 0
with gzip.open(f'{indir}/varcon.SNPs.vcf.gz', 'rt') as f:
    for i, line in enumerate(f):
        if 'CHROM' in line:
            header_line = i
            break
# 1. Load without forcing an index column
vcf = pd.read_csv(f'{indir}/varcon.SNPs.vcf.gz', sep='\t', 
                 compression='gzip', header=header_line, index_col=0)

vcf.columns = [c.replace('#', '').strip() for c in vcf.columns]

def extract_dp(info_string):
    if pd.isna(info_string): return 0
    for item in str(info_string).split(';'):
        if item.startswith('DP='):
            return int(item.split('=')[1])
    return 0

vcf['DP'] = vcf['INFO'].apply(extract_dp)
dp = vcf['DP'].to_numpy()

consistent = [mmread(f'{indir}/{donor}.consistent.mtx.gz').tocsr()
              for donor in donors]

n_barcodes = barcodes.shape[0]


with open(f'{pkldir}/c1_dict.pkl', 'rb') as f:
    c1_dict = pickle.load(f)

with open(f'{pkldir}/c2_dict.pkl', 'rb') as f:
    c2_dict = pickle.load(f)

with open(f'{pkldir}/i1_dict.pkl', 'rb') as f:
    i1_dict = pickle.load(f)

with open(f'{pkldir}/i2_dict.pkl', 'rb') as f:
    i2_dict = pickle.load(f)

###

bcs_donors = [[f'{bc}_{donor}' for donor in donors] for bc in barcodes]
c1_indices = [[c1_dict[x] for x in y] for y in bcs_donors]
c2_indices = [[c2_dict[x] for x in y] for y in bcs_donors]
i1_indices = [[i1_dict[x] for x in y] for y in bcs_donors]
i2_indices = [[i2_dict[x] for x in y] for y in bcs_donors]

### opportunity to filter variants

n_bcs = barcodes.shape[0]

with ProcessPoolExecutor() as ex:
    results = list(tqdm(ex.map(process_bc, range(n_bcs)),total=n_bcs))

c1_counts, c2_counts, i1_counts, i2_counts = map(list, zip(*results))

c1_df = pd.DataFrame(c1_counts, columns=donors, index=barcodes)
c2_df = pd.DataFrame(c2_counts, columns=donors, index=barcodes)
i1_df = pd.DataFrame(i1_counts, columns=donors, index=barcodes)
i2_df = pd.DataFrame(i2_counts, columns=donors, index=barcodes)

###

cov = str(cov_thresh)
c1_df.to_csv(f'{outdir}/{run_name}_c1_df.csv', sep='\t', header=True, index=True)
c2_df.to_csv(f'{outdir}/{run_name}_c2_df.csv', sep='\t', header=True, index=True)
i1_df.to_csv(f'{outdir}/{run_name}_i1_df.csv', sep='\t', header=True, index=True)
i2_df.to_csv(f'{outdir}/{run_name}_i2_df.csv', sep='\t', header=True, index=True)
