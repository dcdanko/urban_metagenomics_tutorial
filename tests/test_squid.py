
from unittest import TestCase
from capalyzer.packet_parser import NCBITaxaTree
from os.path import join, dirname

from squid.squid import (
    root_to_leaf,
    lowest_common_ancestor,
    classify_seq,
)
from squid.api import load_database, classify_reads

KMER_TABLE = join(dirname(__file__), 'small_annotated_31mer_table.csv')
FASTQ_READS = join(dirname(__file__), 'r1.fq')
TAXA_TREE = NCBITaxaTree.parse_files()


class TestSquid(TestCase):

    def test_load_database(self):
        with open(KMER_TABLE) as db_file:
            database = load_database(db_file)
        self.assertEqual(len(database), 100)
        self.assertIn(
            'Mycobacterium marinum M',
            database['AATACGTCCGGAGTATCGACGCACACATGGT']
        )

    def test_classify_reads(self):
        with open(KMER_TABLE) as db_file:
            database = load_database(db_file)
        reads = list(classify_reads(FASTQ_READS, database, TAXA_TREE))
        self.assertEqual(len(reads), 5)

    def test_root_to_leaf(self):
        in_weights = {
            'Staphylococcus aureus': 1,
            'Escherichia coli': 4,
            'Staphylococcus': 10,
        }
        rtl_weights = root_to_leaf(in_weights, TAXA_TREE)
        self.assertEqual(rtl_weights['Staphylococcus aureus'], 11)
        self.assertEqual(rtl_weights['Escherichia coli'], 4)

    def test_lca(self):
        ancestor = lowest_common_ancestor(
            'Staphylococcus aureus',
            'Staphylococcus lugdunensis',
            TAXA_TREE
        )
        self.assertEqual('Staphylococcus', ancestor)

    def test_classify_seq(self):
        with open(KMER_TABLE) as db_file:
            database = load_database(db_file)
        classification = classify_seq(
            'AATACGTCCGGAGTATCGACGCACACATGGT', database, TAXA_TREE
        )
        self.assertIn('Mycobacterium marinum M', classification)
