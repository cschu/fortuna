#!/usr/bin/env python
import sys
import multiprocessing as mp
import subprocess as sub
from collections import namedtuple

from fortuna import findORFs, reverseComplement
import ktoolu_io

BlastHSP = namedtuple('BlastHSP', 'query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq'.split(' '))
ShortBlastHSP = namedtuple('ShortBlastHSP', 'pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq'.split(' '))
ORFDataShort = namedtuple('ORFDataShort', 'orfid mod start end frame orf_prob nlen plen cds strand'.split(' '))
ORFResult = namedtuple('ORFResult', 'read hsp qstart qend qlen subject pident sstart send slen qseq sseq qmod smod dmod strand orftype'.split(' '))
ORFData = namedtuple('ORFData', 'orfid seq orfindex mod start end frame orf_prob nlen plen cds strand'.split(' '))

ORF_MOD = {0: '.free', 1: '.headless', 2: '.tailless'}

def processFileMP(_in, nthreads, blast_db=None, blast_cmd=None, blaster=None, minlen=200):
    pool = mp.Pool(processes=nthreads)
    kwargs = {'minlen': minlen, 'blast_db': blast_db, 'blaster': blaster, 'blast_cmd': blast_cmd}
    results = [pool.apply_async(processRecord, args=(record, ), kwds=kwargs) for record in _in]

    for query_results in (p.get() for p in results):
        for res in query_results:
            yield res

def runBlast(qid, qseq, blastcmd, blastdb, blaster):
    pr = sub.Popen(blastcmd % (blaster, blastdb), shell=True, stdin=sub.PIPE, stderr=sub.PIPE, stdout=sub.PIPE)
    out, err = pr.communicate(('>query_%s\n%s\n' % (qid, qseq)).encode())
    if out.decode().strip():
       return BlastHSP(*(out.decode().strip().split('\n')[0].split()))
    return None

def blastCheckORF(orf, blast_db, blast_cmd, blaster='blastn', min_pid=75):
    def isFullHit(qstart, qend, qlen, sstart, send, slen):
        return (qstart, qend) == (1, qlen) and sorted((sstart, send)) == [1, slen]

    orfseq = orf.seq[int(orf.frame):][orf.start:orf.end + 3]
    orf_result = None
    blast_hit = runBlast(orf.orfid, orfseq, blast_cmd, blast_db, blaster)
    if blast_hit is not None and float(blast_hit.pident) >= min_pid:
        qstart, qend, sstart, send = map(int, blast_hit[6:10])
        qlen, slen = map(int, blast_hit[12:14])

        orftype = 'partial'
        if orf.mod == '.full' and isFullHit(qstart, qend, qlen, sstart, send, slen):
            orftype = 'full'

        orf_result = ORFResult(orf.orfid, ORFDataShort(*(orf[2:])), qstart, qend, qlen,
                               blast_hit.subject, float(blast_hit.pident), sstart, send, slen,
                               blast_hit.qseq, blast_hit.sseq, len(blast_hit.qseq)%3, len(blast_hit.sseq)%3,
                               (len(blast_hit.qseq) - blast_hit.sseq.count('-'))%3, blast_hit.sstrand, orftype)
    return orf_result

def processRecord(record, blast_db=None, blast_cmd=None, blaster=None, minlen=200):
    _id = record[0].strip('>').strip('@')
    _seq = record[1].upper().replace('U', 'T')
    foundORFs = list()
    for strand in '+-':
        if strand == '-':
            _seq = reverseComplement(_seq)
        for frame in range(3):
            for i, orf in enumerate(findORFs(_seq[frame:], minlen=minlen, frame='%c%i' % (strand, frame + 1), allowN=False)):
                if orf is None:
                    continue
                # ORFCandidate = namedtuple('ORFCandidate', 'start end length frame prob'.split(' '))
                mod = ORF_MOD.get(i, '.full')
                if mod in ('.tailless', '.free'):
                    nlen = orf.end - orf.start + 1
                else:
                    nlen = (orf.end + 3) - (orf.start + 1) + 1
                plen = nlen / 3
                orfid = _id.strip().replace(' ', '_')
                orf_result = ORFData(orfid, _seq, i, mod, orf.start, orf.end, frame, orf.prob, nlen, plen, len(foundORFs) + 1, strand)
                if blaster is not None:
                    assert blast_db is not None and blast_cmd is not None
                    blast_result = blastCheckORF(orf_result, blast_db, blast_cmd, blaster=blaster)
                    if blast_result is not None:
                        foundORFs.append((orf_result, blast_result))
                else:
                    foundORFs.append((orf_result, None))

    return foundORFs

if __name__ == '__main__':
    processFileMP(ktoolu_io.readFastq(sys.argv[1]), 4)
