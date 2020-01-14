import sys
import numpy as np
import pandas as pd
import itertools as it
import math
import random
from scipy import stats

path = '/Genomics/grid/users/yushit/courses/QCB455/HapMap3/'

def pname(name):
    '''Prepend the path to the filename'''
    return path + '/' + name

def popen(name):
    '''Open file in the path'''
    return open(pname(name))


### SNP DATA
def read_snp(file):
    '''Read a snp file into a pandas dataframe'''
    return(pd.read_table(
        file,
        sep='\s+', # columns are separated by whitespace
        # names of the columns
        names=[None, 'chromosome', 'morgans', 'position', 'ref', 'alt'],
        index_col=0
    ))

def read_snp_pop(pop):
    return(read_snp(pname(pop + '.snp')))

def get_chr_range(chromosome):
    '''Returns the range of positions where SNPs for a chromosome are kept'''
    filt = SNPs.query('chromosome=={}'.format(chromosome))
    start = SNPs.index.get_loc(filt.iloc[0].name)
    stop  = SNPs.index.get_loc(filt.iloc[-1].name) + 1
    return(start, stop)

SNPs = read_snp_pop('HapMap3')

### Individual data
def read_ind(file):
    '''Read an ind file into a pandas dataframe'''
    return pd.read_table(
        file,
        sep='\s+',
        names=[None, 'sex', 'pop'],
        index_col = 0
    )

def read_ind_pop(pop):
    return(read_ind(pname(pop + '.ind')))


### Genotype data
def read_geno(file):
    '''Reads a geno file into a numpy matrix'''
    return(np.genfromtxt(
        file,                # the file
        dtype='uint8',       # read the data in as 1-byte integers
        delimiter=1,         # 1-byte width data
        missing_values=9,    # 9 indicates missing data
        usemask=True         # return a masked array
    ))

def read_geno_pop(pop):
    return read_geno(pname(pop + '.geno'))

def read_geno_pop_chr(pop, chromosome):
    '''Reads a slice of a geno file into a numpy matrix'''
    f = popen(pop + '.geno')      # open the file
    (start, stop) = get_chr_range(chromosome)
    s = it.islice(f, start, stop) # slice the file
    return read_geno(s)

### Calculate F_ST
def calc_af(geno):
    '''Calculate the allele frequencies for each row in the geno'''
    return geno.mean(axis=1).filled(-1) / 2

def calc_N(geno):
    '''Count the number of non-missing genotypes for each SNP in the geno'''
    return (~geno.mask).sum(axis=1)

def calc_stats(geno):
    '''Calls calc_af and calc_N and returns the results as a tuple'''
    return calc_af(geno), calc_N(geno)

def calc_FST(stats1, stats2):
    ps      = (stats1[0] + stats2[0]) /  2
    deltas2 = (stats1[0] - stats2[0]) ** 2
    pqs = ps*(1-ps)
    num = deltas2 - pqs/2*(1./stats1[1] + 1./stats2[1])
    den = 2*pqs

    FST1 = (num / den)[den != 0].mean()
    FST2 = num.sum() / den.sum()
    return FST1, FST2

### Calculate r_2 between two SNPs
def calc_r(X):
    '''Calculate the correlation between rows of masked array X'''
    # copy the genotype matrix into two separate matrices offset by one row
    A = X[:-1].copy() # everything but last row
    B = X[1:].copy() # everything but first row

    # set the mask of both to True if either one of the entries are masked
    A.mask = B.mask = (A.mask | B.mask)

    # calculate number of non-missing values for each pair of rows
    N = (~A.mask).sum(axis=1)

    # calculate the sums
    sum_A = A.sum(axis=1).filled(0)
    sum_B = B.sum(axis=1).filled(0)

    sum_A2 = (A**2).sum(axis=1).filled(0)
    sum_B2 = (B**2).sum(axis=1).filled(0)

    sum_AB = (A*B).sum(axis=1).filled(0)
    # note that in the above we fill in missing data with a 0

    # remove if denominator is 0 (or if all data was missing)
    rel_rows = ((np.sqrt(N*sum_A2 - sum_A**2)*np.sqrt(N*sum_B2 - sum_B**2) != 0))

    return ( (N[rel_rows]*sum_AB[rel_rows] - sum_A[rel_rows]*sum_B[rel_rows])
            / np.sqrt(N[rel_rows]*sum_A2[rel_rows] - sum_A[rel_rows]**2)
            / np.sqrt(N[rel_rows]*sum_B2[rel_rows] - sum_B[rel_rows]**2))

### Conduct Armitage Trend Test in the simplest way
def calc_corr(X, Y):
    '''calculate the correlation between rows of masked arrays X and Y'''
    # Annoying-but-needed data setup:
    # Make sure arguments are 2-dimensional
    X = np.ma.atleast_2d(X)
    Y = np.ma.atleast_2d(Y)

    # get the masked arrays as integers
    Xnotmask = ~np.ma.getmaskarray(X) * 1
    Ynotmask = ~np.ma.getmaskarray(Y) * 1

    # now onto the calculations:
    # calculate N for each combination of X and Y
    N = np.dot(Xnotmask, Ynotmask.T)

    # calculate the sums
    sum_X  = np.ma.dot(X,    Ynotmask.T).filled()
    sum_X2 = np.ma.dot(X**2, Ynotmask.T).filled()

    sum_Y  = np.ma.dot(Xnotmask, Y.T).filled()
    sum_Y2 = np.ma.dot(Xnotmask, Y.T**2).filled()

    # calculate the crossed sum
    # note the "Y.T * 1"
    # This is needed because we store geno objects as 1-byte integers
    # The product can be greater than 256. we need to increase the storage to a float
    sum_XY = np.ma.dot(X, Y.T * 1.).filled()

    r = ( (N*sum_XY - sum_X * sum_Y)
         / np.sqrt(N*sum_X2 - sum_X**2)
         / np.sqrt(N*sum_Y2 - sum_Y**2) )

    return r, N

def calc_ATT(*args):
    r, N = calc_corr(*args)
    return N*(r**2)

def normalize_geno(geno, p = None):
    if p is None: p = calc_af(geno)
    return ( (geno-(2*p)[:,np.newaxis])/np.sqrt(2*p*(1-p))[:,np.newaxis]).filled(0)

def read_geno_pop_islice(pop, start, stop, step):
    '''Reads a slice of a geno file into a numpy matrix'''
    f = popen(pop + '.geno') # open the file
    s = it.islice(f, start, stop, step) # slice the file
    return read_geno(s)

# H-E Regression to estimate heritability
def HEReg(gen, phe):
    M = len(gen)
    A = gen.T.dot(gen) / M
    Aprime = np.triu(A, 1)
    h2g = Aprime.dot(phe).dot(phe) / (Aprime**2).sum()
    return(h2g)



def sim_bias(geno, phe):
    bias_list = []
    af_seq = {2: [0.00,0.50,1.00],
              5: [0.00,0.20,0.40,0.60,0.80,1.00],
              10:[0.00,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,0.95,1.00],
              20:[0.00,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,
                  0.70,0.75,0.78,0.80,0.82,0.85,0.88,0.90,0.95,1.00]}
    k_seq = [1,2,5,10,20]
    for k in k_seq:
        if k == 1:
            bias_list.append(1-(HEReg(normalize_geno(geno), phe)))
        else:
            h2g_temp = []
            for i in range(k):
                snps_af = (calc_af(geno) >= af_seq[k][i]) & (calc_af(geno) < af_seq[k][i+1])
                new_geno_masks = np.ma.mask_or(np.tile(np.transpose([snps_af]),len(geno[1])), geno.mask)
                new_geno = np.ma.array(geno.data, mask=new_geno_masks)
                new_geno_norm = normalize_geno(new_geno)
                h2g_temp.append(HEReg(new_geno_norm, phe))
            bias_list.append(np.mean(h2g_temp)-1)
    return(bias_list)


def simulate_bias(pop, outdir):
    bias = []
    for i in range(500):
        G = read_geno_pop_islice(pop, 10000+i, 260000+i, 125)
        X = normalize_geno(G)
        B = np.repeat([0,0.1],[1900,100])
        Y = X.T.dot(B)
        bias.append(sim_bias(G,Y))
    bias_df = pd.DataFrame(bias)
    bias_df.to_csv(path_or_buf=outdir,index=False)


if __name__ == '__main__':
    pop = sys.argv[1]
    outdir = sys.argv[2]
    simulate_bias(pop, outdir)
