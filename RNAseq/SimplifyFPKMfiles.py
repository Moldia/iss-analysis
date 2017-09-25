import os

for s in range(25, 29):
    for l in range(4):
        filename = os.path.join(r"E:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\fpkm",
                                "S" + str(s) + "_L00" + str(l+1),
                                "genes.fpkm_tracking")
        with open(filename, "r") as f:
            for line in f:
                line = line.rstrip('\n')
                items = line.split('\t')
                if items[-1] == "OK":
                    geneid = items[3]
                    genename = items[4]
                    fpkm = items[9]



