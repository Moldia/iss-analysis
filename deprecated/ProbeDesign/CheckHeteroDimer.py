import os
import subprocess
from openpyxl import load_workbook


def StartNewBlast():
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess

    if NextProcess < len(probename):
        # print NextProcess
        blastf = probename[NextProcess] + '.fasta'

        if not os.path.isfile(probename[NextProcess] + '_blast.txt'):      # skip if blast is already finished and the file exists

            # local database
            blastn_cline = 'blastn' + ' -query ' + '"' + blastf + '"' + ' -db ProbeDB' + \
                           ' -outfmt 10' + \
                           ' -out ' + '"' + probename[NextProcess] + '_blast.txt" -word_size 7 -strand minus'

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
    while (len(Processes) < 6) and (NextProcess < len(probename)):  # more to do and some spare slots
        StartNewBlast()

# read the probe file
probefile = load_workbook(filename = 'C:/Users/Xiaoyan/OneDrive/worky/Probe_all.xlsx')
probefile = probefile['Sheet1']
probename = probefile.columns[3]
probename = [i.value for i in probename[1:]]

probeseq = probefile.columns[6]
probeseq = [i.value for i in probeseq[1:]]


# prepare BLAST database file
os.mkdir('Temp')
os.chdir('Temp')
with open("ProbeAll.fasta", 'w') as f:
    for c, i in enumerate(probename):
        f.write(">%s\n%s\n" % (i, probeseq[c]))
        with open(i + ".fasta", 'w') as fi:
            fi.write(">%s\n%s\n" % (i, probeseq[c]))

os.system("E:/Bioinformatics/blast-2.2.30+/bin/makeblastdb -in ProbeAll.fasta -out ProbeDB -dbtype nucl")


# run first blast
Processes = []
NextProcess = 0
CheckRunningBlast()  # start the max processes running
while len(Processes) > 0:   # still going on
    CheckRunningBlast()

# read blast output and run the second blast
with open("../BlastOverview.csv", 'w') as f:
    f.write("query,hit,%identity,alignment length,mismatches,gap,q.start,q.end,s.start,s.end,evalue,bitscore\n")
    for i in probename:
        fblast = i + "_blast.txt"
        if os.stat(fblast).st_size == 0:
            os.remove(i + ".fasta")
            os.remove(fblast)
        else:
            with open(fblast) as fb:
                for line in fb:
                    f.write(line)
                os.system("blastn -query " + i + ".fasta" + " -db ProbeDB -outfmt 1" + \
                          " -out " + fblast + " -word_size 7 -strand minus")


