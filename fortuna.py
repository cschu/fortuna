#!/usr/bin/env python
'''
Script to find open reading frames in nucleic acid sequences 
provided in Fasta format.
'''
import re
import sys
import itertools as it

import numpy as np

# import analyse_blast as anabl

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2013-2014, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'cschu@darkjade.net'

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
    compl = dict(zip(alphabet, alphabet[::-1]))
    return ''.join([compl[base] 
                    for base in seq.upper().replace('U', 'T')])[::-1]

def findORFs(seqFrame, start='ATG', stop=OCHRE_AMBER_OPAL, minlen=200, frame=1):
    """
    Searches open reading frames in the input sequence frame.    
    Parameter frame is just a label.
    """
    # First, break down NA-sequence into codons
    codons = [seqFrame[i:i+3] for i in xrange(0, len(seqFrame), 3)]
    starts, stops = [], []
    for i, codon in enumerate(codons):
        if codon == start:
            starts.append(i)
        elif codon in stop:
            stops.append(i)    
    # Find all potential full ORFs(uninterrupted (start, stop) combinations).
    # These represent potential full-length transcripts/peptides.
    # ORF-format: (start, end, length[aa|codons], frame)
    fullORFs = sorted([pair + (pair[1] - pair[0] + 1, frame)
                       for pair in it.product(starts, stops) if pair[0] < pair[1]])
    # Now look for the longest unterminated ORF (free ORF)
    # i.e., the first start after the last detected stop
    stops = [-1] + stops
    ORFstarts = [start for start in starts if start > stops[-1]]
    freeORF = (None, None, 0, frame)
    if ORFstarts:
        lengthFreeORF = len(codons) - ORFstarts[0]
        if (lengthFreeORF * 3) >= minlen:
            freeORF = (ORFstarts[0], len(codons) - 1, lengthFreeORF, frame)
        pass
    
    # Check the compatibility of potential full ORFs
    # (i, j) : (i, j + n) => (i, j) survives
    # (i, j) : ((i + n) < j, j) => (i, j) survives
    validORFs = []
    i = 0
    while True:
        if not fullORFs: break 
        activeORF = fullORFs.pop(0)
        validORFs.append(activeORF)
        invalid = []
        for j in xrange(0, len(fullORFs)):
            if fullORFs[j][0] == activeORF[0]:
                invalid.append(j)
        for p in invalid[::-1]:
            fullORFs.pop(p)

    return [ORF for ORF in validORFs if (ORF[2] * 3) >= minlen], freeORF

def main(argv):
    contigs = anabl_getContigsFromFASTA(argv[0])
    minlen = 200
    orffile = open(re.sub('.fa(sta)?$', '.ORF%i.fa' % minlen, argv[0]), 'wb')
    no_orffile = open(re.sub('.fa(sta)?$', '.NOORF%i.fa' % minlen, argv[0]), 'wb')
    orfseqs = open(re.sub('.fa(sta)?$', '.ORFSEQ%i.fa' % minlen, argv[0]), 'wb')
    
    medianFullORF, medianLongestORF = [], []
    validORFs, validFullORFs = 0.0, 0.0
    i = 0
    for id_, seq in contigs:
        i += 1
        rcSeq = reverseComplement(seq)
        orf_id = id_.strip('>').strip()
        
        # scan all 6 frames
        ORFs, fullORFs, freeORFs = [], [], []
        for frame in xrange(3):
            currentORFs = findORFs(seq[frame:], minlen=minlen, frame=frame + 1)
            fullORFs.extend(currentORFs[0])
            freeORFs.append(currentORFs[1])
            currentORFs = findORFs(rcSeq[frame:], minlen=minlen, frame=frame + 1)
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
            orffile.write('>%s_%if\n%s\n' % (orf_id, len(fullORFs), seq.strip()))
            orfseqs.write('>%s\n%s\n' % (orf_id + 'fORF=%i,%i,frame=%i' % (longestFullORF[0],
                                                                           longestFullORF[1],
                                                                           longestFullORF[3]), seq))
            sys.stdout.write('\r%i sequences processed, %i sequences with ORF, %i sequences with fullORF' % (i, validORFs+validFullORFs, validFullORFs))

        elif longestFreeORF[2] > 0 and longestFreeORF[2] > longestFullORF[2]:
            medianLongestORF.append(longestFreeORF[2] * 3)
            validORFs += 1.0
            orffile.write('>%s_%i\n%s\n' % (orf_id, longestFreeORF[2], seq.strip()))            
            orfseqs.write('>%s\n%s\n' % (orf_id + 'ORF=%i,frame=%i' % (longestFreeORF[0], 
                                                                       longestFreeORF[3]), seq))
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
    minlen = 200

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
