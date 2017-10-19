with open(r"E:\RNAseq\probeseq\bam_tophat_GRCm38\accepted_hits_reformat.bed", 'w') as fw:
    with open(r"E:\RNAseq\probeseq\bam_tophat_GRCm38\accepted_hits.bed", 'r') as f:
        for line in f:
            inline = line.split()
            for item in inline:
                fw.write("%s\t" % item)
            fw.write('\n')
