
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import numpy as np
from sklearn import decomposition
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

plt.close('all')

TT = pd.read_csv('A:\\InSituSequencing\\160417\\Cell_blobs_raw.csv').iloc[1:,1:-2]
T = TT.drop(['NNNN', 'Slc1a2'],1)

counts=T.sum(1)
good_cells = T[counts>=5]

## Make dendrogram for genes:
plt.figure('Genes')
plt.clf()
gene_mat = good_cells.as_matrix().T
distance = ssd.pdist(gene_mat, 'cosine')
linx_genes = sch.linkage(distance, method='weighted')
sch.dendrogram(linx_genes, labels=good_cells.columns, orientation='left', leaf_font_size=12, color_threshold=.9)

plt.figure('Correlations')
plt.clf()
leaves_genes = sch.leaves_list(linx_genes)
corr_mat = good_cells.corr().as_matrix()
plt.imshow(corr_mat[leaves_genes,:][:,leaves_genes], interpolation='none', cmap='hot')
plt.xticks(range(good_cells.shape[1]), good_cells.columns[leaves_genes],rotation='vertical',verticalalignment='top')
plt.yticks(range(good_cells.shape[1]), good_cells.columns[leaves_genes])
plt.grid('on', color=[.1,.1,.1])
plt.colorbar()
plt.clim(0,.2)


## now for cells:
cf = plt.figure('Cells')
cf.clf()
cell_mat = good_cells.as_matrix()
distance_cells = ssd.pdist(cell_mat, 'cosine')
linx_cells = sch.linkage(distance_cells, method='weighted')
#sch.dendrogram(linx, orientation='left', leaf_font_size=12, color_threshold=.9, truncate_mode='level', p=10)

# do PCA
clusters_cells = sch.fcluster(linx_cells,.84, criterion='distance')
pca = decomposition.PCA(n_components=3)
pca.fit(cell_mat)
y = pca.transform(cell_mat)
ax = cf.add_subplot(111, projection='3d')
ax.scatter(y[:,0], y[:,1],y[:,2],c=clusters_cells, cmap='Accent')

## now see which genes are in which clusters
cluster_means = good_cells.groupby(clusters_cells).aggregate(np.mean)

plt.figure('Clusters')
plt.clf()
plt.imshow(cluster_means.as_matrix().T[leaves_genes,:], interpolation='none', cmap='hot')
plt.yticks(range(cluster_means.shape[1]), cluster_means.columns[leaves_genes])
plt.grid('on', color=[.1,.1,.1])
plt.colorbar()
plt.clim(0,2)
plt.xlabel('Cluster')
plt.title('Mean expression')