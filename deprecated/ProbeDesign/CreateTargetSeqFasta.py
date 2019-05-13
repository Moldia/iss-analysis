with open(r"E:\RNAseq\probeseq\probe_targets.fasta", "w") as fw:
    with open(r"E:\RNAseq\probeseq\probe_targets.csv", "r") as f:
        for line in f:
            fw.write('>')
            line = line.rstrip('\n')
            for item in line.split(','):
                fw.write('%s\n' % item)

