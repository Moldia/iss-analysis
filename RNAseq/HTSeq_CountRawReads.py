import HTSeq
# import itertools
import collections


gtf_file = HTSeq.GFF_Reader("/home/ronja/SeqDatabases/mm10/gtf/Mus_musculus.GRCm38.90.gtf")

# all exons, genome location and gene ID
exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
# dictionary of gene ID - transcript names
shortnames = {}

for feature in gtf_file:
    if feature.type == "exon":

        # NB! Ensembl gtf only does not have 'chr' in chromosome name, while UCSC (used in alignment) genome files have
        feature.iv.chrom = 'chr' + feature.iv.chrom

        exons[feature.iv] += feature.attr["gene_id"]
        shortnames[feature.attr["gene_id"]] = feature.attr["transcript_name"]


# sanity check: Sst (ENSMUSG00000004366)
print exons[HTSeq.GenomicPosition("chr16", 23890744, "-")]


# count reads mapped to exons
almnt_file = HTSeq.SAM_Reader("/home/ronja/XQ/S25_L001.sorted.sam")
counts = collections.Counter()

for almnt in almnt_file:
    if not almnt.aligned:
        counts["_unmapped"] += 1
        continue

    gene_ids = set()
    for iv, val in exons[almnt.iv].steps():
        gene_ids |= val

    if len(gene_ids) == 1:
        gene_id = list(gene_ids)[0]
        counts[gene_id] += 1
    elif len(gene_ids) == 0:
        counts["_no_feature"] += 1
    else:
        counts["_ambiguous"] += 1

shortnames["_unmapped"] = '_unmapped'
shortnames["_no_feature"] = '_no_feature'
shortnames["_ambiguous"] = '_ambiguous'


# write raw counts to file
with open("/home/ronja/XQ/S25_L001_HTSeq.csv", 'w') as f:
    for gene_id in counts:
        f.write("%s,%s,%d\n" % (gene_id, shortnames[gene_id], counts[gene_id]))

