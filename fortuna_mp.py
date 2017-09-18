#!/usr/bin/env python
import sys
import multiprocessing as mp
import subprocess as sub

from fortuna import findORFs, reverseComplement
import ktoolu_io

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
       return out.decode().strip().split('\n') # [0].split('\t')
    return list()

def blastCheckORF(orf, blast_db, blast_cmd, blaster='blastn'):
    _seq, _mod = orf[1], orf[3]
    start, end, frame = map(int, orf[4:7])
    orfseq = _seq[frame:][start:end + 3]
    blast_hit = runBlast(orf[0], orfseq, blast_cmd, blast_db, blaster)
    orf_result = None
    if blast_hit and float(blast_hit[0].split()[2]) > 75:
        blast_hit = blast_hit[0].split()
        qstart, qend, sstart, send = map(int, blast_hit[6:10])
        qlen, slen = map(int, blast_hit[12:14])

        orftype = 'partial'
        if _mod == '.full' and (qstart, qend, sstart, send) == (1, qlen, 1, slen):
            orftype = 'full'

        orf_result = (orf[0], orf[2:], qstart, qend, qlen,
                      blast_hit[1], blast_hit[2], sstart, send, slen,
                      blast_hit[20], blast_hit[21], len(blast_hit[20])%3, len(blast_hit[21])%3,
                      (len(blast_hit[20]) - blast_hit[21].count('-'))%3, orftype)
    return orf_result


def processORFs(orfdata):
    orf_count, blast_count = 0, 0
    fullORFs, truncORFs = list(), list()
    covered_regions = set()
    for orf in orfdata:
        orf_count += 1
        _seq, mod = orf[1], orf[3]
        start, end, frame = map(int, orf[4:7])
        orfseq = _seq[frame:][start:end + 3]
        # aaseq = translate(orfseq)
        # test.fq, blastn: 86 98 0.8775510204081632
        # test.fq, blastp: 24 98 0.24489795918367346
        # blast_hit = runBlast(_id, aaseq, BLAST_CMD, BLAST_DB, 'blastp')
        blast_hit = runBlast(orf[0], orfseq, BLAST_CMD, BLAST_DB_CDS, 'blastn')
        if not blast_hit or float(blast_hit[0].split()[2]) < 75:
            continue
        # print(blast_hit)
        blast_hit = blast_hit[0].split()
        blast_count += 1

        """['query_m160210_080137_42165_c100957252550000001823219307011657_s1_p0/287/ccs',
            'GLIADIN.OMEGA.0213.Triticum_aestivum.8099',
            '98.592', '71', '0', '1', '189', '258', '1185', '1115', '4.75e-29', '124',
            '258', '1185', '70', '1', '98.59', '1/1', 'N/A', 'N/A',
            'GAATTCCACAAATGATCTAGGGTACACATGCAACTGTGTCCTTGAGTATGACAACA-CCCTTATCTTCTAG',
            'GAATTCCACAAATGATCTAGGGTACACATGCAACTGTGTCCTTGAGTATGACAACAACCCTTATCTTCTAG']
        """
        qstart, qend, sstart, send = map(int, blast_hit[6:10])
        qlen, slen = map(int, blast_hit[12:14])

        orf_result = (orf[0], orf[2:], qstart, qend, qlen,
                      blast_hit[1], blast_hit[2], sstart, send, slen,
                      blast_hit[20], blast_hit[21], len(blast_hit[20])%3, len(blast_hit[21])%3,
                      (len(blast_hit[20]) - blast_hit[21].count('-'))%3)

        if mod == '.full':
            if qstart == 1 and qend == qlen:
                if send == slen and sstart == 1:
                    # this is a full hit
                    covered_regions.add((start, end))
                    fullORFs.append(orf_result)
                    continue

        truncORFs.append(orf_result)

    return fullORFs, truncORFs




def processRecord(record, blast_db=None, blast_cmd=None, blaster=None, minlen=200):
    _id = record[0].strip('>').strip('@')
    _seq = record[1].upper().replace('U', 'T')
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
                orf_result = (_id.strip().replace(' ', '_'), _seq, i, mod, start, end, frame, orf[4], nlen, plen, len(foundORFs) + 1, strand)
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
