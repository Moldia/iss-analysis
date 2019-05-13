
file = 'E:/PROOOJECTS/12_Neuron_mapping/Probe_design/mRNA_sequences_mouse/ntng1.txt'
fileout = 'E:/PROOOJECTS/12_Neuron_mapping/Probe_design/mRNA_sequences_mouse/ntng1_consensus.txt'

headers = []
sequences = []
seq = []
with open(file, 'r') as f:
    for line in f:
        if line[0] == '>':
            headers.append(line.rstrip('\n'))
            sequences.append(seq)
            seq = str()
        else:
            seq += line.rstrip('\n')
    sequences.append(seq)

start = [1]
end = [2182]
with open(fileout, 'w') as f:
    for i in range(len(start)):
        f.write('>Ntng1_consensus_' + str(i+1) + '\n')
        temp = sequences[1][start[i]-1:end[i]]
        while len(temp)>=70:
            f.write('%s\n' % temp[0:70])
            temp = temp[70:]
        else:
            f.write('%s\n\n' % temp)
