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

adata = ad.read_h5ad("../raw/scATAC-Healthy-Hematopoiesis-191120-adata")
umap25 = np.genfromtxt('../processed/Granja_ATAC_umap25.csv', delimiter=',')

# select reference cells based on 25D UMAP of ATAC
g=0.5

S = pairwise_kernels(umap25, metric='rbf',gamma=g)
median = np.median(S, axis=0)
start = int(np.argmin(median))

dd_rcells = from_knn_dist(X=umap25,
                         n_ret=100,
                         start=start,
                         metric='rbf',
                         metric_params={"gamma": g})

S_to_r = pairwise_similarities(umap25,query=dd_rcells,metric='rbf',metric_params={"gamma": g})


# calculate scores from ATAC
escores = calculate_escores(adata.layers['b'], query=dd_rcells, S=S_to_r, ncpu=16)

escores.to_hdf('../processed/Granja_ATAC_escore.h5', key='escores')

# calculate permuted enrichment scores
P = permute(adata.layers['b'])

pscores = calculate_escores(P, query=dd_rcells, S=S_to_r, ncpu=16)

pscores.to_hdf('../processed/Granja_ATAC_pscore.h5', key='pscores')


