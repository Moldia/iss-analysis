
#-----reverse complement-----
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

# import sequences
name = []
seq = []
primer = []
primer_5 = []
with open('C:/Worky/Projects/Neuromics/Probe_design/Targets.csv', 'r') as f:
    for line in f:
        temp = line.strip('\n').split(',')
        name.append(temp[0])
        seq.append(temp[2])
        primer.append(reverse_complement(temp[2][21:]))
        primer_5.append(reverse_complement(temp[2][25:]))

name = name[1:]
primer = primer[1:]
primer_5 = primer_5[1:]


# write to file
with open('C:/Worky/Projects/Neuromics/Probe_design/Primers.csv', 'w') as f:
    for c, i in enumerate(name):
        f.write("Pr_%s_0,%s\n" % (i, primer[c]))
    for c, i in enumerate(name):
        f.write("Pr_%s_5,%s\n" % (i, primer_5[c]))