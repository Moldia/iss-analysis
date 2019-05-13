
from openpyxl import load_workbook


def reverse_complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases


def change_to_rna(seq, rnapos):
    for c, pos in enumerate(rnapos):
        if pos < 0:
            rnapos[c] = len(seq) + pos
    rnapos = sorted(list(set(rnapos)))

    posoverseq = next(c for x in enumerate(rnapos) if x>=len(seq))
    if posoverseq:
        rnapos = rnapos[:posoverseq]
    elif rnapos[0] < len(seq):
        pass

    for pos in rnapos[::-1]:
        seq = seq[:pos] + 'r' + seq[pos:]

    seq = seq.replace('rT', 'rU')

    return seq


wb = load_workbook(filename='E:\Probes\DirectRNA\DirectRNA.xlsx', read_only=True)
ws = wb['cDNAprobes']

with open('E:\Probes\DirectRNA\DirectRNA_RNAtarget.csv', 'w') as f:
    for row in ws:
        if row[0].value is not None:
            if row[0].value == 'ID':
                for cell in row:
                    f.write('%s,' % cell.value)
                f.write('RNAtarget_arm_sequence,RNA_arm_5,RNA_arm_3,RNA_PLP,probe_construct\n')
            else:
                for j, cell in enumerate(row):
                    if j == 11:
                       f.write('%s,' % cell.value.replace(',', ' '))
                       armlength = cell.value.split(' ')
                       totallength = armlength[1][1:-1]
                       armlength = armlength[0].split(',')
                    else:
                        f.write('%s,' % cell.value)
                        if j == 8:
                            plp = cell.value
                        elif j == 9:
                            rc = reverse_complement(cell.value)
                newplp = ''.join((rc[int(armlength[0]):],
                                  plp[int(armlength[0]):int(totallength)-int(armlength[-1])],
                                  rc[:int(armlength[0])]))
                f.write('%s,%s,%s,%s,' % (rc, rc[int(armlength[0]):], rc[:int(armlength[0])], change_to_rna(newplp, [-1])))
                f.write('%s ' % armlength[-1])
                for c in range(len(armlength)-2):
                    f.write('%s ' % armlength[c+1])
                f.write('%s ' % armlength[0])
                f.write('(%s)\n' % totallength)
        else:
            break




