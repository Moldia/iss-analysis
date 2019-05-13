import subprocess
import os.path
import numpy as np
import math


# FUNCTIONS
def ChopSeq(seq, window, step):
    """ Moving window to chop target sequence """
    chopped_seq = []
    while len(seq) >= window:
        chopped_seq.append(seq[0:window])
        seq = seq[step:]
    return chopped_seq

def CalculateTm(seq):
    """ Calculate Tm of a target candidate, nearest neighbor model"""
    NN_list = ChopSeq(seq, 2, 1)
    NN_table = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    NN_endtable = ['A', 'C', 'G', 'T']
    NN_count = np.zeros(16)
    NN_end = np.zeros(4)
    for p, NN in enumerate(NN_table):
        NN_count[p] = NN_list.count(NN)
    end_list = seq[0]
    for p, NN in enumerate(NN_endtable):
        NN_end[p] = end_list.count(NN)
    # numbers below from Sugimoto et al. NAR (1996)
    NN_H = np.array([-8.0, -9.4, -6.6, -5.6, -8.2, -10.9, -11.8, -6.6, -8.8, -10.5, -10.9, -9.4, -6.6, -8.8, -8.2, -8.0])
    NN_S = np.array([-21.9, -25.5, -16.4, -15.2, -21.0, -28.4, -29.0, -16.4, -23.5, -26.4, -28.4, -25.5, -18.4, -23.5, -21.0, -21.9])
    NN_end_H = np.array([.6, .6, .6, .6])
    NN_end_S = np.array([-9.0, -9.0, -9.0, -9.0])

    sum_H = np.sum(np.multiply(NN_count, NN_H)) + np.sum(np.multiply(NN_end, NN_end_H))
    sum_S = np.sum(np.multiply(NN_count, NN_S)) + np.sum(np.multiply(NN_end, NN_end_S))
    Tm = (sum_H * 1000)/(sum_S + (1.9872 * math.log(1e-7))) - 273.15       # oligo concentration: 1e-7 M
    sum_salt = 0.075 + (3.795 * 0.01**0.5)      # monovalent: 0.075 M, bivalent: 0.01 M
    Tm += 16.6 * math.log10(sum_salt)     # salt correction
    Tm -= 0.72 * 20       # formamide correction
    return Tm

def StartNewBlast():
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess

    if NextProcess < len(sequence):
        blastf = 'Temp/' + sequence[NextProcess][0] + '.fasta'

        # local database
        blastn_cline = 'blastn' + ' -query ' + '"' + blastf + '"' + ' -db E:/Bioinformatics/blast-2.2.30+/db/mouse.all.rna' + \
                       ' -outfmt 10' + \
                       ' -out ' + '"' + 'Temp/' + sequence[NextProcess][0] + '_blast.txt" -word_size 7 -strand plus'

        blastn_process = subprocess.Popen(blastn_cline, shell=True)

        NextProcess += 1
        Processes.append(blastn_process)

def CheckRunningBlast():
    """ Check any running processes and start new ones if there are spare slots """
    global Processes
    global NextProcess

    for p in range(len(Processes)-1,-1,-1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while (len(Processes) < 6) and (NextProcess < len(sequence)):  # more to do and some spare slots
        StartNewBlast()

def ReadBlast(file):
    """ Read the results from blast """
    scores = np.empty((0,10), dtype='float64')
    hits = []
    with open(file, 'r') as fh:
        for l in fh:
            columns = l.split('|')
            if columns[3][0:2] == 'NM' or columns[3][0:2] == 'NR':     # check if the hit is transcript (NM) or non-coding RNA (NR), remove all predicted (XM)
                hits.append(columns[3])
                score = np.array(columns[-1][1:-1].split(','))
                score = score.astype(np.float)
                scores = np.vstack((scores, score))
    uni, idx = np.unique(scores[:,8], return_index=True)
    hits = np.array(hits)
    return (hits[idx], scores[idx,:])

#=======================================================================================================================
# read and format sequences
sequence = []
with open('Book1.csv') as f:
    for line in f:
        line = line.rstrip('\n')
        line = line.split(',')
        sequence.append([line[0],line[1]])

# write query sequences to separate files
os.mkdir('Temp')
Tm = []
with open('mRNA_targets.fasta', 'w') as f:
    for i in range(len(sequence)):
        f.write(">%s\n%s\n\n" %(sequence[i][0], sequence[i][1]))
        Tm.append(CalculateTm(sequence[i][1]))
        with open('Temp/' + sequence[i][0] + '.fasta', 'w') as f2:
            f2.write(">%s\n" % sequence[i][0])
            f2.write("%s\n" % sequence[i][1])

# run blast
Processes = []
NextProcess = 0
CheckRunningBlast()
while len(Processes) > 0:
    CheckRunningBlast()
print("Blast finished!")

# read blast
with open('BlastFirst.csv', 'w') as f1, open('BlastSecond.csv', 'w') as f2:
    f1.write("query,hit,identities,hit_start,hit_end,e_value\n")
    f2.write("query,hit,identities,query_start,query_end,e_value\n")
    blast = [np.zeros((len(sequence),5)), np.zeros((len(sequence),5))]
    for i in range(len(sequence)):
        hits, scores = ReadBlast('Temp/' + sequence[i][0] + '_blast.txt')
        if len(scores):
            if scores[0,0]==100 and scores[0,1]==len(sequence[i][1]):
                f1.write("%s,%s,%d/%d,%d,%d,%s\n"
                        % (sequence[i][0],hits[0],scores[0,1]-scores[0,2]-scores[0,3],scores[0,1],scores[0,6],scores[0,7],repr(scores[0,8])))
                scores = np.delete(scores,0,0)
                hits = np.delete(hits,0,0)
            else:
                f1.write("%s,NA,NA,nan,nan,nan\n" % sequence[i][0])
            if len(scores):
                f2.write("%s,%s,%d/%d,%d,%d,%s\n"
                        % (sequence[i][0],hits[0],scores[0,1]-scores[0,2]-scores[0,3],scores[0,1],scores[0,4],scores[0,5],repr(scores[0,8])))
            else:
                f2.write("%s,NA,NA,nan,nan,nan\n" % sequence[i][0])
        else:
            f1.write("%s,NA,NA,nan,nan,nan\n" % sequence[i][0])
            f2.write("%s,NA,NA,nan,nan,nan\n" % sequence[i][0])


# write Tm values
with open('Tm.csv', 'w') as f:
    for i, T in enumerate(Tm):
        f.write("%s,%s\n" % (sequence[i][0], repr(T)))
