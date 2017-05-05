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

    # the freeORF is a potential coding sequence missing both start and stop codon
    # this can only occur if there are neither starts nor stops present in the sequence
    freeORF = None
    if not starts and not stops:
        freeORF = (0, len(seqFrame), len(seqFrame), frame, 1.0)
    yield freeORF


    # Extract the headless ORF in the sequence,
    # i.e., the sequence from the beginning of the sequence until the first stop.
    # This ORF only exists if it does not contain an AUG, otherwise
    # it would overlap the first full ORF.
    headlessORF = None
    # starts = [Codon()] + starts
    starts = starts + [Codon()]
    stops = stops + [Codon(pos=starts[0].pos + 1)]
    if starts[0].pos > stops[0].pos and stops[0].pos > minlen:
        headlessORF = (0, stops[0].pos, stops[0].pos, frame, 1.0)
        pass
    yield headlessORF
    # Now look for the longest unterminated ORF (taillessORF)
    # i.e., the first start after the last detected stop
    # starts = starts[1:]
    starts = starts[:-1]
    stops = [Codon()] + stops[:-1]
    ORFstarts = (start for start in starts if start.pos > stops[-1].pos)
    taillessORF = None # (-1, -1, 0, 1, 1.0)
    taillessORFStart = None
    try:
        taillessORFStart = next(ORFstarts)
    except:
        pass
    if taillessORFStart is not None:
        lengthTaillessORF = len(seqFrame) - taillessORFStart.pos#n_codons - freeORFStart.pos
        if lengthTaillessORF  >= minlen:
            taillessORF = (taillessORFStart.pos, len(seqFrame), lengthTaillessORF, frame, taillessORFStart.prob)
        pass
    yield taillessORF


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


def processFile(_in, minlen=200):
    orf_mod = {0: '.free', 1: '.headless', 2: '.tailless'}
    for record in _in:
        _id = record[0].strip('>').strip('@')
        _seq = record[1].upper().replace('U', 'T')
        _CDS = 1
        for strand in '+-':
            if strand == '-':
                _seq = reverseComplement(_seq)
            for frame in range(3):
                for i, orf in enumerate(findORFs(_seq[frame:], minlen=minlen, frame='%c%i' % (strand, frame + 1), allowN=False)):
                    if not orf:
                        continue

                    mod = orf_mod.get(i, '')
                    start, end = orf[0], orf[1]
                    nlen = (end + 3) - (start + 1) + 1
                    plen = nlen / 3
                                        
                    yield _id, i, mod, start + 1, end + 3, orf[3], orf[4], nlen, plen, _CDS
                    _CDS += 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--minlen', type=int, default=200)
    ap.add_argument('seqfile', type=str)
    args = ap.parse_args()

    assert(os.path.exists(args.seqfile))
    assert(args.minlen > 0)

    _file, minlen = args.seqfile, args.minlen
    with open(_file + '.orf%i.fa' % minlen, 'w') as orf_out, open(_file + '.pep%i.fa' % minlen, 'w') as pep_out:
        for orf in processFile(readFasta(_file), minlen=minlen):
            # orf: _id, i, mod, start + 1, end + 3, orf[3], orf[4], nlen, plen, _CDS
            _id = orf[0].strip().replace(' ', '_')
            i, mod, start, end, frame, pr, nlen, plen, _CDS = orf[1:]

            head = '%i%s:%i-%i:%i:%s:%.3f' % (i, mod, start + 1, end + 3,  nlen, frame, pr)
            orf_out.write('>%s_CDS%i:%s\n%s\n' % (_id, _CDS, head, _seq[frame:][start:end + 3]))
            head = '%i%s:%i-%i:%i:%s:%.3f' % (i, mod, start + 1, end + 3,  plen, frame, pr)
            pep_out.write('>%s_CDS%i:%s\n%s\n' % (_id, _CDS, head, translate(_seq[frame:][start:end + 3])))
