with open(r"E:\RNAseq\probeseq\probe_targets\accepted_hits.USCS.bed", 'w') as fw:
    with open(r"E:\RNAseq\probeseq\probe_targets\accepted_hits.bed", 'r') as f:
        for line in f:
            fw.write("chr" + line)

