import os
from openpyxl import load_workbook


# all genes that we have probes against
allgenes = []

wb = load_workbook(filename=r"C:\Users\qxyyx\OneDrive\worky\allProbes.xlsm", read_only=True)
ws = wb['cDNAprobes']

for row in ws.rows:
    allgenes.append(row[4].value)

allgenes = list(set(allgenes))
allgenes.remove("name")
# if reference gtf file is from ENSEMBL
# allgenes[allgenes.index("Lphn2")] = "Adgrl2"
# if reference gtf file is from NCBI
allgenes[allgenes.index("Ndnf")] = "A930038C07Rik"
allgenes[allgenes.index("Lamp5")] = "6330527O06Rik"


for s in range(25, 29):
    for l in range(4):
        outfilename = os.path.join(r"E:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\GappedMapping_GRCm38\fpkm_cufflinks",
                                   "S" + str(s) + "_L00" + str(l+1),
                                   "fpkm_sipmle.tsv")

        with open(outfilename, 'w') as fout:
            fout.write("gene_id\tgene_short_name\tfpkm\n")

            filename = os.path.join(r"E:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\GappedMapping_GRCm38\fpkm_cufflinks",
                                    "S" + str(s) + "_L00" + str(l+1),
                                    "genes.fpkm_tracking")
            with open(filename, 'r') as fin:
                for line in fin:
                    line = line.rstrip('\n')
                    items = line.split('\t')
                    if items[-1] == "OK" and items[4] in allgenes:
                        fout.write("%s\t%s\t%f\n" % (items[3], items[4], float(items[9])))

# isoforms file
for s in range(25, 29):
    for l in range(4):
        outfilename = os.path.join(r"E:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\GappedMapping_GRCm38\fpkm_cufflinks",
                                   "S" + str(s) + "_L00" + str(l+1),
                                   "isoforms.fpkm_sipmle.tsv")

        with open(outfilename, 'w') as fout:
            fout.write("track_id\tgene_short_name\tlength\tcoverage\tfpkm\n")

            filename = os.path.join(r"E:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\GappedMapping_GRCm38\fpkm_cufflinks",
                                    "S" + str(s) + "_L00" + str(l+1),
                                    "isoforms.fpkm_tracking")
            with open(filename, 'r') as fin:
                for line in fin:
                    line = line.rstrip('\n')
                    items = line.split('\t')
                    if items[-1] == "OK" and items[4] in allgenes:
                        fout.write("%s\t%s\t%s\t%f\t%f\n" % (items[0], items[3], items[7], float(items[8]), float(items[9])))

