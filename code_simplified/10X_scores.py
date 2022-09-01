# calculate SEMITONES enrichment scores

import warnings
from SEMITONES.enrichment_scoring import calculate_escores
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.enrichment_scoring import permute
from sklearn.preprocessing import binarize
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_kernels

warnings.filterwarnings('ignore')

adata_atac = ad.read_h5ad('../processed/10X_multiome/pbmc_atac.h5ad')

# binarize
adata_atac.layers['b_atac'] = binarize(adata_atac.X)

wnn_umap25= np.genfromtxt('../processed/10X_multiome/wnn_umap25.csv', delimiter=',')


# select reference cells based on 25D UMAP of WNN
g=0.68

S = pairwise_kernels(wnn_umap25, metric='rbf',gamma=g)
median = np.median(S, axis=0)
start = int(np.argmin(median))

dd_rcells = from_knn_dist(X=wnn_umap25,
                         n_ret=100,
                         start=start,
                         metric='rbf',
                         metric_params={"gamma": g})

S_to_r = pairwise_similarities(wnn_umap25,query=dd_rcells,metric='rbf',metric_params={"gamma": g})


# calculate scores from ATAC
n_cpu = 7

escores = calculate_escores(adata_atac.layers['b_atac'], query=dd_rcells, S=S_to_r, ncpu=n_cpu)

escores.to_hdf('../processed/10X_multiome/ATAC_escore.h5', key='escores')

# calculate permuted enrichment scores
P = permute(adata_atac.layers['b_atac'])

pscores = calculate_escores(P, query=dd_rcells, S=S_to_r, ncpu=n_cpu)

pscores.to_hdf('../processed/10X_multiome/ATAC_pscore.h5', key='pscores')


