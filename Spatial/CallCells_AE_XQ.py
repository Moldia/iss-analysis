import numpy as np
from scipy import spatial
import os

input_file = r'G:\Salamander project\New\Cellprofiler\Ch030\Decoding\QT_0.5_details_noNNNN.csv'
output_folder = r'G:\Salamander project\New\Cellprofiler\Analysis'
output_suffix = r'Ch030'

## variables   ##
radius = 30                                             # maximum radius analyzed from an target
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
keepThese = housek + litgen + group2 + group6 + group7
idxKeep = [c for c, gene in enumerate(uNames) if gene in keepThese]
name = name[np.in1d(idxName, idxKeep)]
pos = pos[np.in1d(idxName, idxKeep), :]

# unique names and name index after subsetting
uNames, idxName = np.unique(name, return_inverse=True)

# catagorize genes
typeNames = ['housekeeping', 'literature', 'group2', 'group6', 'group7']
genes = [housek] + [litgen] + [group2] + [group6] + [group7]
geneTypes = [np.tile(c, len(markers)) for c, markers in enumerate(genes)]
genes = [gene for markers in genes for gene in markers]
geneTypes = [c for group in geneTypes for c in group]
uNamesType = [np.where(np.array(genes)==gene)[0][0] for gene in uNames]
uNamesType = np.array(geneTypes)[uNamesType]
idxGroup = uNamesType[idxName]

# build a KDTree to enable fast neighbor search
tree = spatial.KDTree(pos)
# find all neighbors within radius between including all reads
neighbors = tree.query_ball_tree(tree, radius)
AllCells = []
with open(os.path.join(output_folder, 'CellTypeDetails_' + output_suffix + '.csv'), 'w') as f:
    f.write('CellTypeName,Members\n')

    # use only reads from CellType X as query
    for i in range(2, 5):
        idxAllQueryInOneGroup = np.nonzero(idxGroup==i)[0]
        allFoundInOneGroup = []
        querySpot = []
        for query in idxAllQueryInOneGroup:
            NN = np.array(neighbors[query])
            NNGroup = idxGroup[NN]

            # continue only if more than cell_call_limit reads of the same group found
            if np.count_nonzero(NNGroup==i) >= cell_call_limit:
                # calculate distance of a query to each hit
                NNDist = spatial.distance.cdist(np.asmatrix(pos[query,:]), pos[NN,:])[0]

                # sort based on distance
                idx = np.argsort(NNDist)
                NNSortedGroup = NNGroup[idx]
                NNSorted = NN[idx]

                # abort_limit-th closest outsider (invading outsider)
                try:
                    invadingOutsider = \
                        np.nonzero(np.in1d(NNSortedGroup, np.setdiff1d(range(2,5), i)))[0][abort_limit-1]
                except IndexError:      # no outsider within radius
                    invadingOutsider = len(NNSorted)

                # number of insiders closer than the invading outsider
                remainedInsiders = np.nonzero(NNSortedGroup[:invadingOutsider] == i)[0]

                # continue only if more than cell_call_limit neighbors of same group appear before the invading outsider
                if invadingOutsider+1 > cell_call_limit and remainedInsiders.size >= cell_call_limit:
                    # minimum group members (i.e. remove extra fillers farther away than cell_call_limit insiders)
                    minset = min(remainedInsiders[-1], invadingOutsider-1)
                    allFoundInOneGroup.append(tuple(np.append(np.sort(NNSorted[:minset+1]), NNDist[idx][minset])))
                    querySpot.append(NNSorted[0])

        # take unique combinations (INCLUDING HOUSEKEEPING GENS!!!)
        uAllFoundInOneGroup = set(allFoundInOneGroup)
        idx = []
        for group in list(uAllFoundInOneGroup):
            idx.append(allFoundInOneGroup.index(group))

        # organize data and write details to file
        CellsOfOneType = np.zeros((len(uAllFoundInOneGroup), 6+max(geneTypes)+1))
        CellsOfOneType[:,0]= i
        for c, group in enumerate(list(uAllFoundInOneGroup)):
            members = list(group[:-1])
            f.write('%s' % typeNames[i])
            for member in members:
                f.write(',%d' % member)
            f.write('\n')

            members = [a.astype(int) for a in members]
            # query spot XY
            CellsOfOneType[c,1:3] = pos[querySpot[idx[c]],:]
            # recalculated group centroid XY
            CellsOfOneType[c,3:5] = np.mean(pos[members,:], axis=0)
            # max distance
            CellsOfOneType[c,5] = group[-1]
            members = idxGroup[members]

            # different types within a group
            hist = []
            for group in range(0, max(geneTypes)+1):
                hist.append(np.count_nonzero(members==group))
            CellsOfOneType[c,6:] = hist

        AllCells.extend(CellsOfOneType)


# write simplified results
headerline = 'CellTypeName,QueryX,QueryY,CenterX,CenterY,MaxDistance'
for celltype in typeNames:
    headerline = headerline + ',' + celltype
np.savetxt(os.path.join(output_folder, 'CellType_' + output_suffix + '.csv'),
           AllCells, delimiter=',', fmt='%.2f', header=headerline)

