#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script to find open reading frames in nucleic acid sequences
provided in Fasta format.
'''
import re
import sys
import itertools as it


# import analyse_blast as anabl

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2013-2016, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'cschu1981@gmail.com'

OCHRE_AMBER_OPAL = set(['TAA', 'TAG', 'TGA'])

def anabl_getContigsFromFASTA(fn):
    """
    Returns iterator to FASTA file.
    Originates from 'anabl' - BLAST analysing tool, hence the prefix.
    """
    head, seq = None, ''
    for line in open(fn):
        if line[0] == '>':
            if head is not None:
                yield (head, seq)
            head, seq = line.strip().replace(' ', ':').strip('>'), ''
        else:
            seq += line.strip()
    yield (head, seq)

def reverseComplement(seq, alphabet='ACGT'):
    """
    Returns the reverse complement of nucleic acid seqence input.
    """
    # compl = dict(zip(alphabet, alphabet[::-1]))
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]

def calcCodonProbability(codon):
    assert len(codon) == 3
    return 0.25 ** codon.upper().count('N')

class Codon(object):
    def __init__(self, pos, prob):
        self.pos = pos
        self.prob = prob

def findORFs(seqFrame, start='ATG', stop=OCHRE_AMBER_OPAL, minlen=200, frame=1, allowN=True):
    """
    Searches open reading frames in the input sequence frame.
    Parameter frame is just a label.
    """
    start_re = re.compile('[AN][TN][GN]')
    stop_re = re.compile('[TN](([AN][AGN])|([GN][AN]))')
    # First, break down NA-sequence into codons
    codons = ((i, seqFrame[i:i+3]) for i in xrange(0, len(seqFrame), 3))
    starts, stops = list(), list()
    p_start, p_stop = list(), list()
    for i, codon in codons:
        if codon == start or (allowN and start_re.match(codon)):
            #starts.append(i)
            #p_start.append(calcCodonProbability(codon))
            starts.append(Codon(i, calcCodonProbability(codon)))
        elif codon in stop or (allowN and stop_re.match(codon)):
            #stops.append(i)
            #p_stop.append(calcCodonProbability(codon))
            stops.append(Codon(i, calcCodonProbability(codon)))
    n_codons = i + 1
    # Find all potential full ORFs(uninterrupted (start, stop) combinations).
    # These represent potential full-length transcripts/peptides.
    # ORF-format: (start, end, length[aa|codons], frame)
    """
    fullORFs = sorted([pair + (pair[1] - pair[0] + 1, frame)
                       for pair in it.product(starts, stops) if pair[0] < pair[1]])
    """
    fullORFs = sorted((pair[0].pos, pair[1].pos, pair[1].pos - pair[0].pos, frame, pair[0].prob * pair[1].prob)
                      for pair in it.product(starts, stops) if pair[0].pos < pair[1].pos)
    # Now look for the longest unterminated ORF (free ORF)
    # i.e., the first start after the last detected stop
    stops = [-1] + stops
    ORFstarts = (start for start in starts if start.pos > stops[-1].pos)
    freeORF = (-1, -1, -1, -1, -1)
    freeORFStart = None 
    try: 
        freeORFStart = next(ORFstarts) 
    except:
        pass
    # freeORF = (None, None, 0, frame)
    if freeORFStart is not None:
        lengthFreeORF = len(seqFrame) - freeORFStart.pos#n_codons - freeORFStart.pos
        if lengthFreeORF  >= minlen:
            freeORF = (freeORFStart.pos, len(seqFrame), lengthFreeORF, frame, freeORFStart.prob)
        pass
    yield freeORF

    # Check the compatibility of potential full ORFs
    # (i, j) : (i, j + n) => (i, j) survives
    # (i, j) : ((i + n) < j, j) => (i, j) survives
    validORFs = []
    i = 0

    while fullORFs:
        activeORF = fullORFs.pop(0)
        #validORFs.append(activeORF)
        if activeORF[2] >= minlen:
            yield activeORF
        invalid = list()
        for j in xrange(0, len(fullORFs)):
            if fullORFs[j][0] == activeORF[0]:
                invalid.append(j)
            elif fullORFs[j][1] == activeORF[1]:
                invalid.append(j)
            elif fullORFs[j][0] <= activeORF[1]:
                invalid.append(j)
        for p in invalid[::-1]:
            fullORFs.pop(p)

    # return [ORF for ORF in validORFs if (ORF[2] * 3) >= minlen], freeORF

def main(argv):
    import numpy as np
    contigs = anabl_getContigsFromFASTA(argv[0])
    minlen = 0 #200
    orffile = open(re.sub('.fa(sta)?$', '.ORF%i.fa' % minlen, argv[0]), 'wb')
    no_orffile = open(re.sub('.fa(sta)?$', '.NOORF%i.fa' % minlen, argv[0]), 'wb')
    orfseqs = open(re.sub('.fa(sta)?$', '.ORFSEQ%i.fa' % minlen, argv[0]), 'wb')

    medianFullORF, medianLongestORF = [], []
    validORFs, validFullORFs = 0.0, 0.0
    i = 0
    for id_, seq in contigs:
        if 'N' in seq:
            continue
        i += 1
        rcSeq = reverseComplement(seq)
        orf_id = id_.strip('>').strip()

        # scan all 6 frames
        ORFs, fullORFs, freeORFs = [], [], []
        for frame in xrange(3):
            currentORFs = findORFs(seq[frame:], minlen=minlen, frame=frame + 1)
            fullORFs.extend(currentORFs[0])
            freeORFs.append(currentORFs[1])
            currentORFs = findORFs(rcSeq[frame:], minlen=minlen, frame=-(frame + 1))
            fullORFs.extend(currentORFs[0])
            freeORFs.append(currentORFs[1])

        # get longest ORFs
        # ORF-format: (start, end, length[aa|codons], frame)
        longestFullORF = (None, None, 0, None)
        longestFreeORF = (None, None, 0, None)
        if fullORFs:
            longestFullORF = sorted(fullORFs, key=lambda x:x[2])[-1]
        if freeORFs:
            longestFreeORF = sorted(freeORFs, key=lambda x:x[2])[-1]

        if longestFullORF[2] > 0 and longestFullORF[2] >= longestFreeORF[2]:
            medianFullORF.append(longestFullORF[2] * 3)
            validFullORFs += 1.0
            if longestFullORF[3] > 0:
                lFullSeq = seq[longestFullORF[3] - 1:][longestFullORF[0] * 3:longestFullORF[1] * 3 + 1]
            else:
                lFullSeq = rcSeq[abs(longestFullORF[3]) - 1:][longestFullORF[0] * 3:longestFullORF[1] * 3 + 1]
            orffile.write('>%s_%if\n%s\n' % (orf_id, len(fullORFs), seq.strip()))
            orfseqs.write('>%s\n%s\n' % (orf_id + 'fORF=%i,%i,frame=%i' % (longestFullORF[0],
                                                                           longestFullORF[1],
                                                                           longestFullORF[3]), lFullSeq))
            sys.stdout.write('\r%i sequences processed, %i sequences with ORF, %i sequences with fullORF' % (i, validORFs+validFullORFs, validFullORFs))

        elif longestFreeORF[2] > 0 and longestFreeORF[2] > longestFullORF[2]:
            medianLongestORF.append(longestFreeORF[2] * 3)
            validORFs += 1.0
            if longestFreeORF[3] > 0:
                lFreeSeq = seq[longestFreeORF[3] - 1:][longestFreeORF[0] * 3:-1]
            else:
                lFreeSeq = rcSeq[abs(longestFreeORF[3]) - 1:][longestFreeORF[0] * 3:-1]
            orffile.write('>%s_%i\n%s\n' % (orf_id, longestFreeORF[2], seq.strip()))
            orfseqs.write('>%s\n%s\n' % (orf_id + 'ORF=%i,frame=%i' % (longestFreeORF[0],
                                                                       longestFreeORF[3]), lFreeSeq))
            sys.stdout.write('\r%i sequences processed, %i sequences with ORF, %i sequences with fullORF' % (i, validORFs+validFullORFs, validFullORFs))

        else:
            no_orffile.write('>%s_%i\n%s\n' % (orf_id, 0, seq.strip()))


            # continue
        #sys.stdout.write('\r%i sequences processed, %i sequences with ORF (%.3f, medLen=%f), %i sequences with fullORF (%.3f, medLen=%f (%i))' % (i, validORFs+validFullORFs, (validORFs+validFullORFs)/i,
        #                                                                                                                                         np.median(medianFullORF+medianLongestORF),
        #                                                                                                                                         validFullORFs, validFullORFs/i, np.median(medianFullORF), len(medianFullORF)))
        #break
        pass

    sys.stdout.write('\r%i sequences processed, %i sequences with ORF (%.3f, medLen=%f), %i sequences with fullORF (%.3f, medLen=%f (%i))' % (i, validORFs+validFullORFs, (validORFs+validFullORFs)/i, np.median(medianFullORF+medianLongestORF), validFullORFs, validFullORFs/i, np.median(medianFullORF), len(medianFullORF)))



    orfseqs.close()
    orffile.close()
    no_orffile.close()
    sys.stdout.write('\n')




def main2(argv):

    # print findORFs('ATGATGACGGGAAAACCACCCGCGTGATGACGTGA')
    contigs = anabl_getContigsFromFASTA(argv[0])
    minlen = 0#200

    orffile = open(re.sub('.fa(sta)?$', '.ORF%i.fa' % minlen, argv[0]), 'wb')
    no_orffile = open(re.sub('.fa(sta)?$', '.NOORF%i.fa' % minlen, argv[0]), 'wb')
    orfseqs = open(re.sub('.fa(sta)?$', '.ORFSEQ%i.fa' % minlen, argv[0]), 'wb')


    i, validORFs, validFullORFs = 0, 0, 0
    medianFullORF, medianLongestORF = [], []
    for id_, seq in contigs:
        i += 1
        print i, id_, seq
        rcSeq = reverseComplement(seq)
        orfFound = False

        ORFs = []
        fullORFs, longestORF= [], []
        for frame in xrange(3):
            ORFs = findORFs(seq[frame:], minlen=minlen, frame=frame+1)
            fullORFs.extend(ORFs[0])
            longestORF.append((ORFs[1][1], frame+1))
            ORFs = findORFs(rcSeq[frame:], minlen=minlen, frame=-(frame+1))
            fullORFs.extend(ORFs[0])
            longestORF.append((ORFs[1][1], -(frame+1)))


        longestFullORF, fullFrame = sorted(fullORFs, key=lambda x:(x[1]-x[0]+1))[-1] if len(fullORFs) > 0 else (None, None)

        #longestFullORF, fullFrame = sorted([((ORF[1]-ORF[0]+1), ORF[2])
        #                                    for ORF in fullORFs])[-1] if len(fullORFs) > 0 else (0, None)

        # longestORF = (ORFstarts[0], len(codons) - ORFstarts[0])
        if longestORF:
            longestORF, frame_ = sorted(longestORF, key=lambda x:x[1])[-1]
            lengthLongestORF = longestORF
        else:
            lengthLongestORF = 0

        if longestFullORF is None:
            lengthLongestFullORF = 0
        else:
            lengthLongestFullORF = longestFullORF[1] - longestFullORF[0] + 1

        if lengthLongestFullORF > 0 and lengthLongestFullORF >= lengthLongestORF:
            medianFullORF.append(lengthLongestFullORF * 3)
            validFullORFs += 1.0
            orffile.write('>%s_%if\n%s\n' % (id_.strip('>').strip(),
                                             sum(map(len, fullORFs)),
                                             seq.strip()))
            orf_id = id_.strip('>').strip() + 'fORF=%i,%i,frame=%i' % longestFullORF + (fullFrame,)
            orfseqs.write('>%s\n%s\n' % (orf_id, seq))
            sys.stdout.write('\r%i sequences processed, %i sequences with ORF, %i sequences with fullORF' % (i, validORFs+validFullORFs, validFullORFs))

        elif lengthLongestORF > 0 and lengthLongestORF > lengthLongestFullORF:
            medianLongestORF.append(lengthLongestORF * 3)
            validORFs += 1.0
            orffile.write('>%s_%i\n%s\n' % (id_.strip('>').strip(),
                                            longestORF[0], seq.strip()))
            sys.stdout.write('\r%i sequences processed, %i sequences with ORF, %i sequences with fullORF' % (i, validORFs+validFullORFs, validFullORFs))
        else:
            no_orffile.write('>%s_%i\n%s\n' % (id_.strip('>').strip(),
                                               0, seq.strip()))
            orf_id = id_.strip('>').strip() + 'ORF=%i,frame=%i' % (longestORF, frame_)
            orfseqs.write('>%s\n%s\n' % (orf_id, seq))
            # continue
        # break
        pass

    print 'X'
    #sys.stdout.write('%i sequences processed, %i sequences with ORF (%.3f, medLen=%f), %i sequences with fullORF (%.3f, medLen=%f (%i))' % (i, validORFs+validFullORFs, (validORFs+validFullORFs)/i,
    #                                                                                                                                        np.median(medianFullORF+medianLongestORF),
    #                                                                                                                                        validFullORFs, validFullORFs/i, np.median(medianFullORF), len(medianFullORF)))

    orfseqs.close()
    orffile.close()
    no_orffile.close()
    sys.stdout.write('\n')



    pass


if __name__ == '__main__': main(sys.argv[1:])

"""AA,CTA,GAT,TTT,TCT,TAT,AAA,TAG
ATG,CTT,CTT,ATT,CAT,TCA,ACG,TCA,CCA,AAC,ATG,GTA,TCA,TTT,CTC,CTA,CTT,TGT,TAT,CCA,AAT,TTG,TAT,CTA,AAG,AAA,TTT,ATT,GTT,TAC,ACA,TTT,ATC,TAA,TAA,CGT,CTA,AAC,ATG,GTC,AAC,GGG,GCA,TTA,ATT,GGG,TGA,TTAACGAAGATGTAGCCCTTTGTAAGGCGTGGGTAATTGTTAACGAAGACAATGTCAATGGAATGTACCAAGCATCATAATATTTTGGGGCTTGAGTTTATGCTTATTTCAACAGCAACCTACAAATTTTGTCAGGCAAAGAAACCCACATAGGGGGAGGCATTGAGTCTCGCTTGAAA
"""
