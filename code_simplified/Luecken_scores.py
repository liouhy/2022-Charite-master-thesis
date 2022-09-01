# calculate SEMITONES enrichment scores

import warnings
from SEMITONES.enrichment_scoring import calculate_escores
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.enrichment_scoring import permute
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_kernels

warnings.filterwarnings('ignore')

adata = ad.read_h5ad("../processed/Luecken_multiome_BMMC-r_adata.h5ad")

wnn_umap25= np.genfromtxt('../processed/Luecken_wnn_umap25.csv', delimiter=',')

iATAC = np.where(adata.var['feature_types']=='ATAC')[0]

# select reference cells based on 25D UMAP of WNN
g=0.6

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

escores = calculate_escores(adata.X[:,iATAC], query=dd_rcells, S=S_to_r, ncpu=n_cpu)

escores.to_hdf('../processed/Luecken_ATAC_escore.h5', key='escores')

# calculate permuted enrichment scores
P = permute(adata.X[:,iATAC])

pscores = calculate_escores(P, query=dd_rcells, S=S_to_r, ncpu=n_cpu)

pscores.to_hdf('../processed/Luecken_ATAC_pscore.h5', key='pscores')


