import numpy as np
from scipy import spatial


input_file = r'C:\Users\Xiaoyan\SkyDrive\worktemp\160408_160222UCL_GCaMP-cortex2\Benchmarking\CP_161018\Decoding\QT_0.35_0.0001_details.csv'

## variables   ##
radius = 100                                             # maximum radius analyzed from an target
cell_call_limit = 2                                      # number of self targets to call a cell type
abort_limit = 1

# number of other targets to rule-out cell type
##    lists    ##
housek = ['40S','60S','atp5','eif2']			                # housekeeping genes
group2 = ['DMBT1','MARCS','AMPN','Copia','Gypsy','LINE.2']               # markers for group 2 - expressing TEs (Nuria)
group6 = ['NDUA7','PELO','unkG6','ANGL2']                                # markers for group 6 - expressing ECM genes (Ahmed)
group7 = ['SRGN','CYTF','LYSC','unkG7','CCR4', 'RORC', 'CD11b','MPEG1']  # markers for group 7 - probably immunecells  (Connie)
litgen = ['KAZD','CIRBP','Sox9']                                         # genes targeted based on literature
unwant = ['NNNN']					                # for filtering out (like NNNN)
target_group = group2

# import data
data = np.genfromtxt(input_file, dtype=None, delimiter=',', skip_header=1, usecols=[1,2,3], names=['name', 'posX', 'posY'])
name = data['name']
pos = np.transpose(np.vstack((data['posX'], data['posY'])))
del data

# remove reads not needed for neighbor search
uNames, idxName = np.unique(name, return_inverse=True)
# remove unwanted genes
# idxUnwanted = [c for c, gene in enumerate(uNames) if gene in unwant]
# name = name[np.in1d(idxName, idxUnwanted, invert=True)]
# pos = pos[np.in1d(idxName, idxUnwanted, invert=True),:]

# or keep only the ones needed in neighbor search
keepThese = housek + group2 + group6 + group7 + litgen
idxKeep = [c for c, gene in enumerate(uNames) if gene in keepThese]
name = name[np.in1d(idxName, idxKeep)]
pos = pos[np.in1d(idxName, idxKeep), :]

# unique names and name index
uNames, idxName = np.unique(name, return_inverse=True)

# catagorize genes
genes = [housek] + [litgen] + [group2] + [group6] + [group7]
geneGroup = [np.tile(c, len(markers)) for c, markers in enumerate(genes)]
genes = [gene for markers in genes for gene in markers]
geneGroup = [c for group in geneGroup for c in group]
uNamesGroup = [np.nonzero(gene==genes)[0] for gene in uNames]
idxGroup = uNamesGroup[idxName]

# KDTree search
dist = spatial.distance_matrix(pos, pos)
tree = spatial.KDTree(pos)
neighbors = tree.query_ball_tree(tree, radius)
Cells = []
for i in range(2, 5):
    queryAllInOneGroup = np.nonzero(idxGroup==i)
    for query in queryAllInOneGroup[0]:
        NN = np.array(neighbors[query])
        NNdist = dist[query, NN]
        idx = np.argsort(NNdist)
        NN = NN[idx]




