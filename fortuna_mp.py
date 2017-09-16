#!/usr/bin/env python
import sys
import multiprocessing as mp

from fortuna import findORFs, reverseComplement
import ktoolu_io

ORF_MOD = {0: '.free', 1: '.headless', 2: '.tailless'}

def processFileMP(_in, nthreads, minlen=200):
    pool = mp.Pool(processes=nthreads)

    for query_results in [pool.apply(processRecord, args=(record, ), kwds={'minlen': minlen}) for record in _in]:
        for res in query_results:
            yield res
    

def processRecord(record, minlen=200):
    _id = record[0].strip('>').strip('@')
    _seq = record[1].upper().replace('U', 'T')
    _CDS = 1
    foundORFs = list()
    for strand in '+-':
        if strand == '-':
            _seq = reverseComplement(_seq)
        for frame in range(3):
            for i, orf in enumerate(findORFs(_seq[frame:], minlen=minlen, frame='%c%i' % (strand, frame + 1), allowN=False)):
                if not orf:
                    continue

                mod = ORF_MOD.get(i, '.full')
                start, end = orf[0], orf[1]
                if mod in ('.tailless', '.free'):
                    nlen = end - start + 1
                else:
                    nlen = (end + 3) - (start + 1) + 1
                plen = nlen / 3
                foundORFs.append((_id.strip().replace(' ', '_'), _seq, i, mod, start, end, frame, orf[4], nlen, plen, _CDS, strand))
                _CDS += 1
    return foundORFs

processFileMP(ktoolu_io.readFastq(sys.argv[1]), 4)
