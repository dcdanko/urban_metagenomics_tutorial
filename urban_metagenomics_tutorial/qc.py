
from gimmebio.kmers import make_kmers
from math import log2


def kmer_entropy(seqs, k=31):
    tbl = {}
    for seq in seqs:
        for kmer in make_kmers(seq, k, canon=True):
            tbl[kmer] = 1 + tbl.get(kmer, 0)
    H, total = 0, sum(tbl.values())
    for count in tbl.values():
        p = count / total
        H += p * log2(p)
    return -H


def gc_content(seqs):
    gc, at = 0, 0
    for seq in seqs:
        for base in seq:
            if base in ['G', 'C', 'g', 'c']:
                gc += 1
            else:
                at += 1
    return gc / (gc + at)
