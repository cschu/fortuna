#!/usr/bin/env python

import re
import sys
import itertools as it

import numpy as np

import analyse_blast as anabl



OCHRE_AMBER_OPAL = set(['TAA', 'TAG', 'TGA'])


def reverseComplement(seq, alphabet='ACGT'):
    compl = dict(zip(alphabet, alphabet))
    return ''.join([compl[base] for base in seq.upper().replace('U', 'T')])

def findORFs(seq, start='ATG', stop=OCHRE_AMBER_OPAL, minlen=200, frame=1):
    codons = [seq[i:i+3] for i in xrange(0, len(seq), 3)]
    starts, stops = [], []
    for i, codon in enumerate(codons):
        if codon == start:
            starts.append(i)
        elif codon in stop:
            stops.append(i)
    fullORFs = sorted([pair for pair in it.product(starts, stops) if pair[0] < pair[1]])
    if len(stops) > 0:
        ORFstarts = [start for start in starts if start > stops[-1]]
        longestORF = None, 0
        if len(ORFstarts) > 0:
            longestORF = (ORFstarts[0], len(codons) - ORFstarts[0])
        if (longestORF[1] * 3) < minlen:
            longestORF = None, 0
    else:
        longestORF = None, 0
        
    validORFs = []
    i = 0
    while True:
        if len(fullORFs) == 0: break
        activeORF = fullORFs.pop(0)
        validORFs.append(activeORF)
        invalid = []
        for j in xrange(0, len(fullORFs)):
            # if (start, end1) and (start, end2): 
            # remove (start, max(end1,end2)) since min(end1, end2) terminates the ORF
            if fullORFs[j][0] == activeORF[0]:
                invalid.append(j)
        for p in invalid[::-1]:
            fullORFs.pop(p)

    return [ORF + (frame,) 
            for ORF in validORFs 
            if (ORF[1]-ORF[0]+1) * 3 >= minlen], longestORF



def checkFrames(seq, minlen=200):
    revcomp = reverseComplement(seq)
    allORFs = [(1, findORFs(seq, minlen=minlen)), 
               (2, findORFs(seq[1:], minlen=minlen)), 
               (3, findORFs(seq[2:], minlen=minlen)), 
               (-1, findORFs(revcomp, minlen=minlen)), 
               (-2, findORFs(revcomp[1:], minlen=minlen)), 
               (-3, findORFs(revcomp[2:], minlen=minlen))]
    return [ORFs for ORFs in allORFs if len(ORFs[1]) > 0]



def main(argv):
    
    # print findORFs('ATGATGACGGGAAAACCACCCGCGTGATGACGTGA')
    contigs = anabl.getContigsFromFASTA(argv[0])
    minlen = 200

    orffile = open(re.sub('.fa(sta)?$', '.ORF%i.fa' % minlen, argv[0]), 'wb')
    no_orffile = open(re.sub('.fa(sta)?$', '.NOORF%i.fa' % minlen, argv[0]), 'wb')
    orfseqs = open(re.sub('.fa(sta)?$', '.ORFSEQ%i.fa' % minlen, argv[0]), 'wb')


    i, validORFs, validFullORFs = 0, 0, 0
    medianFullORF, medianLongestORF = [], []
    for id_, seq in contigs:
        i += 1
        rcSeq = reverseComplement(seq)
        orfFound = False

        ORFs = []
        fullORFs, longestORF = [], []
        for frame in xrange(3):            
            ORFs = findORFs(seq[frame:], minlen=minlen, frame=frame+1)
            fullORFs.extend(ORFs[0])
            longestORF.append((ORFs[1][1], frame+1))
            ORFs = findORFs(rcSeq[frame:], minlen=minlen, frame=-(frame+1))
            fullORFs.extend(ORFs[0])
            longestORF.append((ORFs[1][1], -(frame+1)))

        
        longestFullORF, fullFrame = sorted(fullORFs, key=lambda x:(x[1]-x[0]+1))[-1] if len(fillORFs) > 0 else (None, None)

        #longestFullORF, fullFrame = sorted([((ORF[1]-ORF[0]+1), ORF[2]) 
        #                                    for ORF in fullORFs])[-1] if len(fullORFs) > 0 else (0, None)

        # longestORF = (ORFstarts[0], len(codons) - ORFstarts[0])
        longestORF, frame_ = sorted(longestORF, key=lambda x:x[1])[-1]

        lengthLongestFullORF = longestFullORF[1] - longestFullORF[0] + 1
        lengthLongestORF = longestORF[1]        

        if lengthLongestFullORF > 0 and lengthLongestFullORF >= lengthLongestORF:
            medianFullORF.append(lengthLongestFullORF * 3)
            validFullORFs += 1.0
            orffile.write('>%s_%if\n%s\n' % (id_.strip('>').strip(), 
                                             sum(map(len, fullORFs)), 
                                             seq.strip()))
            orf_id = id_.strip('>').strip() + 'fORF=%i,%i,frame=%i' % longestFullORF + (fullFrame,)
            orfseqs.write('>%s\n%s\n' % (orf_id, seq))            

        elif lengthLongestORF > 0 and lengthLongestORF > lengthLongestFullORF:
            medianLongestORF.append(lengthLongestORF * 3)
            validORFs += 1.0
            orffile.write('>%s_%i\n%s\n' % (id_.strip('>').strip(), 
                                            longestORF[0], seq.strip()))
        else:
            no_orffile.write('>%s_%i\n%s\n' % (id_.strip('>').strip(), 
                                               0, seq.strip()))
            orf_id = id_.strip('>').strip() + 'ORF=%i,frame=%i' % (longestORF[0], frame_)
            orfseqs.write('>%s\n%s\n' % (orf_id, seq))            
            # continue
        sys.stdout.write('\r%i sequences processed, %i sequences with ORF (%.3f, medLen=%f), %i sequences with fullORF (%.3f, medLen=%f (%i))' % (i, validORFs+validFullORFs, (validORFs+validFullORFs)/i, 
                                                                                                                                                 np.median(medianFullORF+medianLongestORF),
                                                                                                                                                 validFullORFs, validFullORFs/i, np.median(medianFullORF), len(medianFullORF)))       
        pass
    
    
    orfseqs.close()
    orffile.close()
    no_orffile.close()
    sys.stdout.write('\n')



    pass


if __name__ == '__main__': main(sys.argv[1:])       


