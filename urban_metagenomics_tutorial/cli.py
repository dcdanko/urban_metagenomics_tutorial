
import click

from Bio import SeqIO
from capalyzer.packet_builder.sub_factories.toxonomy_long_form import longform_taxa

from .qc import kmer_entropy, gc_content


def subset_fastq(fastq, n):
    seqs = []
    for rec in SeqIO.parse(fastq, 'fastq'):
        seqs.append(rec.seq)
        if len(seqs) >= n:
            break
    return seqs


@click.group()
def main():
    pass


@main.command('qc-stats')
@click.option('-n', '--number-reads', type=int, deafult=10000)
@click.option('-o', '--outfile', type=click.File('w'), default='-')
@click.argument('fastqs', nargs=-1, type=click.File('r'))
def qc_stats(n, outfile, fastqs):
    """Print QC Stats for each file."""
    print('filename\tgc_content\tkmer_entropy', file=outfile)
    for fastq in fastqs:
        seqs = subset_fastq(fastq, n)
        entropy = kmer_entropy(seqs)
        gc = gc_content(seqs)
        print(f'{fastq.name}\t{gc}\t{entropy}', file=outfile)


@main.command('parse-taxa')
@click.option('-o', '--outfile', type=click.File('w'), default='-')
@click.argument('reports', nargs=-1, type=click.File('r'))
def parse_taxa(outfile, reports):
    """Reparse kraken tables into a more useful format."""
    tbl = longform_taxa(reports)
    tbl.to_csv(outfile)
