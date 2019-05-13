__author__ = 'Xiaoyan'
# padlock probe design pipeline
# need to input mRNA sequences
# Tm calculation with nearest neighbor model
# 2016-2-7


import sys
import subprocess
import os.path
import datetime
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
    for c, NN in enumerate(NN_table):
        NN_count[c] = NN_list.count(NN)
    end_list = seq[0]
    for c, NN in enumerate(NN_endtable):
        NN_end[c] = end_list.count(NN)
    # numbers below from Sugimoto et al. NAR (1996)
    NN_H = np.array([-8.0, -9.4, -6.6, -5.6, -8.2, -10.9, -11.8, -6.6, -8.8, -10.5, -10.9, -9.4, -6.6, -8.8, -8.2, -8.0])
    NN_S = np.array([-21.9, -25.5, -16.4, -15.2, -21.0, -28.4, -29.0, -16.4, -23.5, -26.4, -28.4, -25.5, -18.4, -23.5, -21.0, -21.9])
    NN_end_H = np.array([.6, .6, .6, .6])
    NN_end_S = np.array([-9.0, -9.0, -9.0, -9.0])

    sum_H = np.sum(np.multiply(NN_count, NN_H)) + np.sum(np.multiply(NN_end, NN_end_H))
    sum_S = np.sum(np.multiply(NN_count, NN_S)) + np.sum(np.multiply(NN_end, NN_end_S))
    Tm = (sum_H * 1000)/(sum_S + (1.9872 * math.log(1e-7))) - 273.15       # oligo concentration: 1e-7 M
    sum_salt = 0.075 + (3.795 * 0.01**0.5)      # monovalent: 0.075 M, bivalent: 0.01 M
    Tm += 16.6 * math.log10(sum_salt)       # salt correction
    Tm -= 0.72 * 20       # formamide correction
    return Tm


def StartNewBlast():
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess

    if NextProcess < len(sites):
        # print NextProcess
        blastf = fname + '_' + str(NextProcess+1) + '.fasta'

        # local database
        blastn_cline = 'blastn' + ' -query ' + blastf + ' -db E:/Bioinformatics/blast-2.2.30+/db/refseq_rna.00' + \
                       ' -outfmt 10' + \
                       ' -out ' + fname + '_' + str(NextProcess+1) + '_blast.txt'

        # # remotely run blastn
        # blastn_cline = 'blastn' + ' -query ' + blastf + ' -db refseq_rna' + \
        #                ' -outfmt 10' + \
        #                ' -out ' + fname + '_' + str(NextProcess+1) + '_blast_new.txt' + ' -remote'

        blastn_process = subprocess.Popen(blastn_cline, shell=True)
        NextProcess += 1
        Processes.append(blastn_process)


def CheckRunning():
    """ Check any running processes and start new ones if there are spare slots """
    global Processes
    global NextProcess

    for p in range(len(Processes)-1,-1,-1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while (len(Processes) < 6) and (NextProcess < len(sites)):  # More to do and some spare slots
        StartNewBlast()


def ReadBlast(file):
    """ Read the results from blast """
    specific = 1
    with open(file, 'r') as fh:
        for l in fh:
            if specific:
                scores = l.split('|')[-1]
                scores = scores.split(',')
                if scores[1] != '100.00':
                    specific = 0
    return specific


def LigationSite(seq, armlength, maxdev):
    """ Determin ligation site """
    left = max(armlength-maxdev, len(seq)-armlength-maxdev)
    right = min(len(seq)-armlength+maxdev, armlength+maxdev)
    Ligation_table = []
    sites = range(left, right+1)
    for i in sites:
        Ligation_table.append([seq[0:i], seq[i:]])
    return Ligation_table


# ======================================================================================================================
# MAIN
print("This is a computational padlock probe design tool!")

# input from user keyboard
fastafile = raw_input("Fasta file containing mRNA sequence(s). Fasta header must be present! : ")
armlength = input("Predefined arm length: ")
if (not isinstance(armlength, int)) or armlength < 0:
    print('Specified armlength is not integer or negative. Terminate.')
    sys.exit()
maxtolerance = input("Maximum deviation of arm length from predefined: ")
if (not isinstance(maxtolerance, int)) or maxtolerance < 0:
    print('Specified armlength is not integer or negative. Terminate.')
    sys.exit()

# read sequence file
try:
    with open(str(fastafile), 'r') as f:
        sequence = [line.rstrip('\n') for line in f]
except IOError:
    print('Input file not found. Terminate.')
    sys.exit()

# find the fasta header lines starting with '>'
headerlines = [i for i in range(len(sequence)) if (sequence[i] and sequence[i][0] == '>')]

# format sequences
sequences = []
headers = []
if len(headerlines) == 0:
    print('No sequence in fasta format found. Terminate.')
    sys.exit()
elif len(headerlines) > 1:
    for i in range(len(headerlines)-1):
        sequences.append(''.join(sequence[headerlines[i]+1:headerlines[i+1]-1]))
        headers.append(sequence[headerlines[i]])

# the last entry (the only one in the case of one sequence in the whole file)
sequences.append(''.join(sequence[headerlines[-1]+1:len(sequence)-1]))
headers.append(sequence[headerlines[-1]])
del sequence


# process sequence by sequence, get target candidates that fulfill Tm requirements, write to file
Tm_table = []
Chopped_site = []
t = str(datetime.datetime.now())[:19]
t = t.translate(None, '-')
t = t.translate(None, ':')
t = t.translate(None, ' ')
os.mkdir('TempFolder' + t)
with open('Tm_selected' + t + '.csv', 'w') as f:
    for i, seq in enumerate(sequences):
        fname = headers[i].split('|')
        fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'     # access number and name
        fname = fname.translate(None, ' ')      # remove whitespace

        f.write("%s\n" % headers[i])
        Tm_temp = []
        Chopped_site_temp = []
        Chopped_length_temp = []

        c = 0
        for j in range(-maxtolerance*2, maxtolerance*2+1):
            chopped_list = ChopSeq(seq, armlength*2+j, 1)
            for k, chopped_sequence in enumerate(chopped_list):
                Tm = CalculateTm(chopped_sequence)
                if 65 < Tm < 75:        # Tm thresholding
                    Tm_temp.append(Tm)
                    f.write("%s,%f,%d,%d\n" % (chopped_sequence, Tm, k+1, k+armlength*2+j))
                    c += 1
                    Chopped_site_temp.append(k)
                    Chopped_length_temp.append(armlength*2+j)

                    # write files that can be used as input in blastn
                    with open(fname + '_' + str(c) +'.fasta', 'w') as fblast:
                        fblast.write(">target_%d\n%s\n" % (c, chopped_sequence))
        Tm_table.append(Tm_temp)
        Chopped_site.append([Chopped_site_temp,Chopped_length_temp])
# del sequences

# blast targets
print ("Starting blast..")
for i, sites in enumerate(Chopped_site):
    sites = sites[0]
    fname = headers[i].split('|')
    fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'
    fname = fname.translate(None, ' ')      # read files written in the previous section

    Processes = []
    NextProcess = 0
    CheckRunning()  # start the max processes running
    while len(Processes) > 0:   # still going on
        # time.sleep(0.1)
        CheckRunning()
print ("Blast finished!")


# read blast files
for i, sites in enumerate(Chopped_site):
    fname = headers[i].split('|')
    fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'
    fname = fname.translate(None, ' ')      # read files written in the previous section

    blast_bw = []
    Tm_temp = np.array(Tm_table[i])
    Chopped_site_temp = np.array(sites)
    for j, target in enumerate(sites[0]):
        fblast = fname + '_' + str(j+1) +'_blast.txt'
        blast_bw.append(ReadBlast(fblast))

    # find sequences that are specific enough
    idx_blast = np.nonzero(blast_bw)[0]
    Chopped_site[i] = Chopped_site_temp[:,idx_blast]
    Tm_table[i] = Tm_temp[idx_blast]


# determine ligation sites, get candidates that fulfill delta Tm requirements
Arm_table = []
Tm_arms = []
for i, sites in enumerate(Chopped_site):
    Arm_temp = []
    Tm_temp = []
    idx_temp = []
    for j in range(len(sites[0,:])):
        seq = sequences[i][sites[0,j]:sites[0,j]+sites[1,j]]
        padlockarms = LigationSite(seq, armlength, maxtolerance)
        for arms in padlockarms:
            Tm_1 = CalculateTm(arms[0])
            Tm_2 = CalculateTm(arms[1])
            if abs(Tm_1 - Tm_2) < 5:        # delta Tm thresholding
                idx_temp.append(j)
                Arm_temp.append(arms)
                Tm_temp.append([Tm_1, Tm_2])
    Arm_table.append(Arm_temp)
    Tm_arms.append(Tm_temp)
    Tm_table[i] = Tm_table[i][idx_temp]
    Chopped_site[i] = sites[:,idx_temp]

# assemble padlock sequences, and write to file
linker_up = "CTGATTCCTTTGACTCACATTGCGTCTATTTAGTGGAGCC"
linker_down = "CTATCTTCTTT"
with open('probes_output' + t + '.csv', 'w') as f:
    f.write("target, Tm, target_5, Tm_5, target_3, Tm_3, start_pos, end_pos, padlock\n")
    for i, armset in enumerate(Arm_table):
        f.write("%s\n" % headers[i])
        for j, arms in enumerate(armset):
            padlock = arms[1] + linker_up + 'XXXX' + linker_down + arms[0]
            f.write("%s,%f,%s,%f,%s,%f,%d,%d,%s\n"
                    % (arms[0]+arms[1], Tm_table[i][j],
                       arms[1], Tm_arms[i][j][1],
                       arms[0], Tm_arms[i][j][0],
                       Chopped_site[i][0,j]+1, Chopped_site[i][0,j]+Chopped_site[i][1,j],
                       padlock))



# TODO: ranking of probes (self-complementary? GC%?)