import os
from openpyxl import load_workbook


# all genes that we have probes against
allgenes = []

wb = load_workbook(filename=r"C:\Users\Xiaoyan\OneDrive\worky\allProbes.xlsm", read_only=True)
ws = wb['Sheet1']

for row in ws.rows:
    allgenes.append(row[4].value)

allgenes = list(set(allgenes))
allgenes.remove("name")


for s in range(25, 29):
    for l in range(4):
        outfilename = os.path.join(r"G:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\UngappedMapping\fpkm",
                                   "S" + str(s) + "_L00" + str(l+1),
                                   "fpkm_sipmle.csv")

        with open(outfilename, 'w') as fout:

            filename = os.path.join(r"G:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\UngappedMapping\fpkm",
                                    "S" + str(s) + "_L00" + str(l+1),
                                    "genes.fpkm_tracking")
            with open(filename, 'r') as fin:
                for line in fin:
                    line = line.rstrip('\n')
                    items = line.split('\t')
                    if items[-1] == "OK" and items[4] in allgenes:
                        fout.write("%s,%s,%f\n" % (items[3], items[4], float(items[9])))


