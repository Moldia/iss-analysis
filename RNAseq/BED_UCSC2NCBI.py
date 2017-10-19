with open(r"E:\RNAseq\probeseq\probe_targets\accepted_hits.mm10.bed", 'w') as fw:
    with open(r"E:\RNAseq\probeseq\probe_targets\accepted_hits.USCS.mm10.bed", 'r') as f:
        for line in f:
            fw.write(line[3:])

