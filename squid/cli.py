import click
from capalyzer.packet_parser import NCBITaxaTree

from .api import load_database, classify_reads


@click.group()
def main():
    pass


@main.command('classify')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('database', type=click.File('r'))
@click.argument('fastq_file', type=click.File('r'))
def cli_classify(outfile, database, fastq_file):
    taxa_tree = NCBITaxaTree.parse_files()
    database = load_database(database)
    classification_generator = classify_reads(fastq_file, database, taxa_tree)
    for read_id, classification in classification_generator:
        print(read_id + ',' + classification, file=outfile)
