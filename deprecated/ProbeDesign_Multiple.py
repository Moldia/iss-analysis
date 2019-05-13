__author__ = 'Xiaoyan'
# padlock probe design pipeline
# need to input mRNA sequences
# Tm calculation with nearest neighbor model
# adapted to multiple probes per target
# 2016-3-24


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


# def LigationSite(seq, armlength, maxdev):
#     """ Determin ligation site """
#     early = max(armlength-maxdev, len(seq)-armlength-maxdev)
#     late = min(len(seq)-armlength+maxdev, armlength+maxdev)
#     sites = range(early, late+1)
#     Ligation_table = []
#     for i in sites:
#         Ligation_table.append([seq[0:i], seq[i:]])
#     return Ligation_table


def StartNewBlast():
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess

    if NextProcess < len(sites):
        # print NextProcess
        blastf = fname + '_' + str(NextProcess+1) + '.fasta'

        if not os.path.isfile(fname + '_' + str(NextProcess+1) + '_blast.txt'):      # skip if blast is already finished and the file exists

            # local database
            blastn_cline = 'blastn' + ' -query ' + '"' + blastf + '"' + ' -db E:/Bioinformatics/blast-2.2.30+/db/mouse.all.rna' + \
                           ' -outfmt 10' + \
                           ' -out ' + '"' + fname + '_' + str(NextProcess+1) + '_blast.txt" -word_size 7 -strand plus'

            # # remotely run blastn
            # blastn_cline = 'blastn' + ' -query ' + blastf + ' -db refseq_rna' + \
            #                ' -outfmt 10' + \
            #                ' -out ' + fname + '_' + str(NextProcess+1) + '_blast_new.txt' + ' -remote'

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
    while (len(Processes) < 6) and (NextProcess < len(sites)):  # more to do and some spare slots
        StartNewBlast()


def ReadBlast(file):
    """ Read the results from blast """
    global funmap
    global Not_mapped
    specific = True
    mappable = False
    with open(file, 'r') as fh:
        for l in fh:
            if specific:
                columns = l.split('|')
                if columns[3][0:2] == 'NM' or columns[3][0:2] == 'NR':     # check if the hit is transcript (NM) or non-coding RNA (NR), remove all predicted (XM)
                    scores = columns[-1].split(',')
                    if 2*armlength*.5 < int(scores[2]) < 2*armlength and float(scores[1]) > 80 and int(scores[5]) < armlength-4 and int(scores[6] > armlength+5):
                        # more than 50% coverage, 80% homology, and non-target sequence covers ligation site +- 5
                        specific = False
                    if not mappable:
                        if float(scores[1])==100 and int(scores[2])==2*armlength:
                            mappable = True
        if not mappable:
            with open(file[0:-10] + '.fasta') as f:
                seq = f.readlines()
                funmap.write('Could not map sequence in ' + file[:-10] + '!\n')
                funmap.write(seq[1])
            Not_mapped.append(int(file[:-10].split('_')[-1])-1)
    return specific


# ======================================================================================================================
# MAIN

# input from user keyboard
fastafile = raw_input("Fasta file containing mRNA sequence(s). Fasta header must be present! : ")
armlength = 20
maxtolerance = 0
interval = 20

# read sequence file
try:
    with open(str(fastafile), 'r') as f:
        sequence = [line.rstrip('\n') for line in f]
        sequence.append([])
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
Tm_list = []
Chopped_site = []
t = str(datetime.datetime.now())[:19]
t = t.translate(None, '-')
t = t.translate(None, ':')
t = t.translate(None, ' ')
os.mkdir('TempFolder' + t)
with open('Tm_selected_' + t + '.csv', 'w') as f:
    for i, seq in enumerate(sequences):
        fname = headers[i].translate(None, ':')
        fname = fname.translate(None, '/')
        fname = fname.split('|')
        fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'     # access number and name
        fname = fname.translate(None, ' ')      # remove whitespace

        f.write("%s\n" % headers[i])
        Tm_temp = []
        Chopped_site_temp = []

        c = 0
        chopped_list = ChopSeq(seq, armlength*2, 1)
        for j, chopped_sequence in enumerate(chopped_list):
            Tm = CalculateTm(chopped_sequence)
            Tm_temp.append(Tm)
            if 65 < Tm < 75:        # Tm thresholding
                f.write("%s,%f,%d\n" % (chopped_sequence, Tm, j+1))
                c += 1
                Chopped_site_temp.append(j)

                # write files that can be used as input in blastn
                with open(fname + '_' + str(c) +'.fasta', 'w') as fblast:
                    fblast.write(">target_%d\n%s\n" % (c, chopped_sequence))

        Tm_list.append(Tm_temp)
        Chopped_site.append(Chopped_site_temp)


# blast targets
print ("Starting blast..")
for i, sites in enumerate(Chopped_site):
    fname = headers[i].translate(None, ':')
    fname = headers[i].translate(None, '/')
    fname = fname.split('|')
    fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'
    fname = fname.translate(None, ' ')      # read files written in the previous section

    Processes = []
    NextProcess = 0
    CheckRunningBlast()  # start the max processes running
    while len(Processes) > 0:   # still going on
        # time.sleep(0.1)
        CheckRunningBlast()
print ("Blast finished!")


# read blast files
Candidate_site = []
with open('Unmappable_' + t + '.txt', 'w') as funmap:
    for i, sites in enumerate(Chopped_site):
        funmap.write("%s\n" % headers[i])

        fname = headers[i].translate(None, ':')
        fname = headers[i].translate(None, '/')
        fname = fname.split('|')
        fname = 'TempFolder' + t + '/' + fname[3] + fname[4] + '_target'
        fname = fname.translate(None, ' ')      # read files written in the previous section
        Not_mapped = []

        blast_bw = []
        Chopped_site_temp = np.array(sites)
        for j, target in enumerate(sites):
            fblast = fname + '_' + str(j+1) + '_blast.txt'
            blast_bw.append(ReadBlast(fblast))

        # find sequences that are specific enough
        idx_blast = np.nonzero(blast_bw)[0]
        site_blast = Chopped_site_temp[idx_blast]

        # write unmappable sites
        Not_mapped = Chopped_site_temp[Not_mapped]
        recorded = False
        for c, j in enumerate(Not_mapped):
            if c == 0:
                funmap.write("\nUnmapped sequence starting position(s):\n%d" % (j+1))
                recorded = True
            else:
                if j == temp+1:
                    funmap.write("-")
                    recorded = False
                    if c == len(Not_mapped)-1:
                        funmap.write("%d" % (j+1))
                else:
                    if recorded:
                        funmap.write(",%d" % (j+1))
                    else:
                        funmap.write("%d,%d" % (temp+1, j+1))
                        recorded = True
            temp = j

        # continuous regions
        if len(site_blast):
            idx_pair_start = np.nonzero(site_blast[1:] - site_blast[0:-1] != 1)[0]
            if len(idx_pair_start) == 0 and len(site_blast):    # only one continuous region exists
                idx_pair_start = np.array([0])
                idx_pair_end = np.array([len(site_blast)-1])
            else:
                if idx_blast[0] == 0:
                    idx_pair_start = np.append(idx_blast[0], np.add(idx_pair_start, 1))
                else:
                    idx_pair_start = np.add(idx_pair_start, 1)

                idx_pair_end = np.zeros(len(idx_pair_start), dtype=np.int)
                for j in range(len(idx_pair_start)-1):
                    p = idx_pair_start[j]
                    c = 0
                    while site_blast[p+c+1] == site_blast[p+c]+1:
                        c += 1
                    idx_pair_end[j] = p+c

            site_pair_start = site_blast[idx_pair_start]
            site_pair_end = site_blast[idx_pair_end]
            site_pair_end[-1] = site_blast[-1]

            Candidate_site.append(np.vstack((site_pair_start, site_pair_end)))
        else:
            Candidate_site.append(np.zeros((2, 0)))
        funmap.write("\n\n")


# design as many probes as possible
# starting from the 3' end of the mRNA molecules
# no overlap between probes
# at least xx nt (default: 20) between two probes
Probes = []
Probes_pos = []
for i, sites in enumerate(Candidate_site):
    probes = []
    probes_pos = []
    try:
        current = Candidate_site[i][1][-1]
        sequence = sequences[i]
        while current >= Candidate_site[i][0][0]:
            probes.append(sequence[current:current+armlength*2])
            probes_pos.append(current)
            current -= armlength*2 + interval
            idx_max_start = np.amax(np.append(np.nonzero(Candidate_site[i][0] <= current), -1))
            idx_min_end = np.amin(np.append(np.nonzero(Candidate_site[i][1] >= current), len(sequence)))
            while idx_max_start != idx_min_end and current >= Candidate_site[i][0][0]:
                current -= 1
                idx_max_start = np.amax(np.nonzero(Candidate_site[i][0] <= current))
                idx_min_end = np.amin(np.nonzero(Candidate_site[i][1] >= current))
    except IndexError:      # no usable fragment
        probes = []
        probes_pos = []
    Probes.append(probes)
    Probes_pos.append(probes_pos)


# assemble padlock sequences
# ligation site in the middle by default
linker_up = "LINKERFIRST"
linker_down = "LINKERSECOND"
Padlocks = []
for i, probes in enumerate(Probes):
    padlocks = []
    for probe in probes:
        padlocks.append(probe[armlength:] + linker_up + 'XXXX' + linker_down + probe[0:armlength])
    Padlocks.append(padlocks)


# write to file
with open('probes_output_' + t + '.csv', 'w') as f:
    f.write("target, Tm, start_pos, end_pos, padlock\n")
    for i, header in enumerate(headers):
        f.write("%s\n" % header)
        for j, probe in enumerate(Probes[i]):
            f.write("%s,%f,%d,%d,%s\n"
                    % (probe, Tm_list[i][Probes_pos[i][j]],
                       Probes_pos[i][j]+1, Probes_pos[i][j]+armlength*2,
                       Padlocks[i][j]))
        f.write("\n")


# TODO: ranking of probes (self-complementary? GC%?)
