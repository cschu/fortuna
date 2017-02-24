#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script to find open reading frames in nucleic acid sequences
provided in Fasta format.
'''
import os
import re
import sys
import argparse
import itertools as it

from ktoolu_io import readFasta


try:
    tmp = xrange(1)
except:
    xrange = range

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2013-2017, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'cschu1981@gmail.com'

OCHRE_AMBER_OPAL = set(['TAA', 'TAG', 'TGA'])

codonWheel = { 'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
               'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
               'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
               'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
               'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
               'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
               'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
               'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
               'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
               'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
               'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
               'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
               'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F', }

wobbleWheel = { 'AC': 'T', 'CC': 'P', 'CG': 'R', 'CT': 'L', 'GC': 'A', 'GG': 'G', 'GT': 'V', 'TC': 'S' }

ambiguityWheel = { 'RAY': 'B', 'TGY': 'C', 'GAY': 'D', 'GAR': 'E', 'TTY': 'F', 'CAY': 'H', 'ATH': 'I', 'AAR': 'K', 'TTR': 'L', 'YTR': 'L', 'AAY': 'N', 'CAR': 'Q', 'AGR': 'R', 'MGR': 'R', 'AGY': 'S', 'TAY': 'Y', 'SAR': 'Z' }

def translate(seq, frame=1):
    assert 1 <= frame <= 3
    return ''.join(translateCodon(seq[frame - 1:][i:i+3].upper()) for i in range(0, len(seq) - 3, 3))

def translateCodon(codon):
    assert len(codon) == 3
    aa = codonWheel.get(codon, 'X')
    if aa == 'X':
        aa = wobbleWheel.get(codon[:2], 'X')
        if aa == 'X':
            aa = ambiguityWheel.get(codon, 'X')
    return aa


def reverseComplement(seq, alphabet='ACGT'):
    """
    Returns the reverse complement of nucleic acid seqence input.
    """
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]

def calcCodonProbability(codon):
    assert len(codon) == 3
    return 0.25 ** codon.upper().count('N')

class Codon(object):
    def __init__(self, pos=-1, prob=1.0):
        self.pos = pos
        self.prob = prob

def findORFs(seqFrame, start='ATG', stop=OCHRE_AMBER_OPAL, minlen=200, frame=1, allowN=True):
    """
    Searches open reading frames in frame one of the input sequence.
    Parameter frame is just for the label.
    """
    start_re = re.compile('[AN][TN][GN]')
    stop_re = re.compile('[TN](([AN][AGN])|([GN][AN]))')
    # First, break down NA-sequence into codons
    codons = ((i, seqFrame[i:i+3]) for i in xrange(0, len(seqFrame), 3))
    starts, stops = list(), list()
    p_start, p_stop = list(), list()
    for i, codon in codons:
        if codon == start or (allowN and start_re.match(codon)):
            starts.append(Codon(i, calcCodonProbability(codon)))
        elif codon in stop or (allowN and stop_re.match(codon)):
            stops.append(Codon(i, calcCodonProbability(codon)))
    n_codons = i + 1
    # Find all potential full ORFs(uninterrupted (start, stop) combinations).
    # These represent potential full-length transcripts/peptides.
    # ORF-format: (start, end, length[aa|codons], frame)
    fullORFs = sorted((pair[0].pos, pair[1].pos, pair[1].pos - pair[0].pos, frame, pair[0].prob * pair[1].prob)
                      for pair in it.product(starts, stops) if pair[0].pos < pair[1].pos)
    # Now look for the longest unterminated ORF (free ORF)
    # i.e., the first start after the last detected stop
    stops = [Codon()] + stops
    ORFstarts = (start for start in starts if start.pos > stops[-1].pos)
    freeORF = None # (-1, -1, 0, 1, 1.0)
    freeORFStart = None
    try:
        freeORFStart = next(ORFstarts)
    except:
        pass
    if freeORFStart is not None:
        lengthFreeORF = len(seqFrame) - freeORFStart.pos#n_codons - freeORFStart.pos
        if lengthFreeORF  >= minlen:
            freeORF = (freeORFStart.pos, len(seqFrame), lengthFreeORF, frame, freeORFStart.prob)
        pass
    yield freeORF


    # The ORFlist is sorted so that
    # (i, j) != (i', j'): i <= i' AND j <= j'
    # Check the compatibility of potential full ORFs
    # (i, j) : (i, j + n) => (i, j) survives
    # (i, j) : ((i + n) < j, j) => (i, j) survives
    validORFs = []
    i = 0

    while fullORFs:
        activeORF = fullORFs.pop(0)
        if activeORF[2] >= minlen:
            yield activeORF
        invalid = list()
        for j in xrange(0, len(fullORFs)):
            if fullORFs[j][0] == activeORF[0]:
                # fullORF[j] starts at activeORF, but is longer,
                # thus it is truncated by activeORF's stop codon
                invalid.append(j)
            elif fullORFs[j][1] == activeORF[1]:
                # fullORF[j] and activeORF end at same position,
                # but activeORF is longer than fullORF[j]
                invalid.append(j)
            elif fullORFs[j][0] <= activeORF[1]:
                # fullORF[j] is contained in activeORF
                invalid.append(j)
        for p in invalid[::-1]:
            fullORFs.pop(p)




if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--minlen', type=int, default=200)
    ap.add_argument('seqfile', type=str)
    args = ap.parse_args()

    assert(os.path.exists(args.seqfile))
    assert(args.minlen > 0)

    _file, minlen = args.seqfile, args.minlen
    with open(_file + '.orf.fa', 'w') as orf_out, open(_file + '.pep.fa', 'w') as pep_out:
        for _id, _seq in readFasta(_file):
            _seq = _seq.upper().replace('U', 'T')
            _CDS = 1
            for strand in '+-':
                if strand == '-':
                    _seq = reverseComplement(_seq)
                for frame in range(3):
                    for i, orf in enumerate(findORFs(_seq[frame:], minlen=minlen, frame='%c%i' % (strand, frame + 1), allowN=False)):
                        mod = ''
                        if i == 0 and orf:
                            mod = '.free'
                        elif not orf:
                            continue

                        start, end = orf[0], orf[1]
                        nlen = (end + 3) - (start + 1) + 1
                        plen = nlen / 3
                        _id = _id.strip().replace(' ', '_')
                        head = '_CDS%i:%i%s:%i-%i:%i:%s:%.3f' % (_CDS, i, mod, start + 1, end + 3,  nlen, orf[3], orf[4])
                        orf_out.write('%s:%s\n%s\n' % (_id, head, _seq[frame:][start:end + 3]))
                        head = '_CDS%i:%i%s:%i-%i:%i:%s:%.3f' % (_CDS, i, mod, start + 1, end + 3,  plen, orf[3], orf[4])
                        pep_out.write('%s:%s\n%s\n' % (_id, head, translate(_seq[frame:][start:end + 3])))
                        _CDS += 1
