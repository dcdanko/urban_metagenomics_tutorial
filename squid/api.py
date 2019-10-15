from Bio import SeqIO

from .squid import classify_seq


def load_database(database_file, sep=','):
    """Return a doictionary of kmers -> taxa.

    Parse the database file into a python dictionary.
    The keys in this dictionary will be DNA sequences
    and the values will be the taxonomic annotation
    where this sequence is found.

    The database file is composed of lines with two
    values each: a kmer and its taxonomic annotation. 
    These lines are seperated by a comma.
    """
    out = {}
    for line in database_file:
        tkns = line.strip().split(sep)
        if len(tkns) == 2:
            out[tkns[0]] = tkns[1]
    return out


def classify_reads(seq_file, database, taxa_tree, filetype='fastq', k=31):
    """Produce tuple of read ids and classification for each read in a file."""
    for rec in SeqIO.parse(seq_file, filetype):
        yield rec.id, classify_seq(rec.seq, database, taxa_tree, k=k)
