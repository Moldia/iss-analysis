dirname = 'C:/Worky/Projects/Neuromics/Probe_design'


# load acronym and full  header file
with open(dirname + '/mouse.acronymheaders.txt', 'r') as f:
    Acronyms = [line.rstrip('\n') for line in f]
with open(dirname + '/mouse.selectedheaders.txt', 'r') as f:
    Headers = [line.rstrip('\n') for line in f]

# load target genes
genes = []
with open("C:/Worky/Projects/Neuromics/NewProbes_Nathan/ListOfProbes_with_targets.csv", 'r') as f:
    for line in f:
        line = line.split(',')
        if line[5] != '0':
            genes.append(line[1])
genes = list(set(genes))

# find hits from the list
hits = []
for gene in genes:
    if gene not in Acronyms:
        hits.append([])
    else:
        hit = [c for c, header in enumerate(Acronyms) if header == gene]
        hits.append(hit)

# load database sequences
Seq = []
with open(dirname + '/mouse.selectedseqs.txt', 'r') as f:
    Seq = [line.rstrip('\n') for line in f]

# retrieve sequences
targets = []
for c, hit in enumerate(hits):
    if len(hit) == 0:
        targets.append([])
    elif len(hit) == 1:
        targets.append(Seq[hit[0]])
    else:
        with open(dirname + '/' + genes[c] + '.fasta', 'w') as f:
            for multi in hit:
                f.write('%s\n' % Headers[multi])
                f.write('%s\n\n' % Seq[multi])
