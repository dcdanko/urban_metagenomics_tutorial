
from .constants import NO_TAXA


def lowest_common_ancestor(taxa1, taxa2, taxa_tree):
    """Return a string giving the lowest common ancestor of taxa{1, 2}.

    Compute the most specific taxonomic rank shared by taxa1 and taxa2.
    e.g. if taxa1 is 'Staphylococcus aureus' and taxa2 is 'Staphylococcus
    lugdunensis' this function should return 'Staphylococcus'.

    Use `taxa_tree.ancestors(<taxon>)` to get a list of ancestors for taxon.
    """
    ancestors1 = set(taxa_tree.ancestors(taxa1))
    for ancestor in taxa_tree.ancestors(taxa2):
        if ancestor in ancestors1:
            return ancestor


def root_to_leaf(input_weights, taxa_tree):
    """Return a dictionary mapping taxa to weights.

    Each output weight should be the sum of input weights
    for a taxa and all of its ancestor taxa up to the 'root'
    of the tree.

    Use `taxa_tree.ancestors(<taxon>)` to get a list of ancestors for taxon.
    """
    output_weights = {}
    for taxon, weight in input_weights.items():
        for ancestor in taxa_tree.ancestors(taxon):
            out_weight = output_weights.get(taxon, 0)
            output_weights[taxon] = out_weight + input_weights.get(ancestor, 0)
    return output_weights


def classify_seq(seq, database, taxa_tree, k=31):
    """Return a taxonomic classification for a sequence."""
    num_kmers = len(seq) - k + 1
    raw_weights = {}
    for i in range(num_kmers):
        kmer = seq[i:i + k]
        try:
            taxon = database[kmer]
            raw_weights[taxon] = 1 + raw_weights.get(taxon, 0)
        except KeyError:  # the kmer is not in the database, common
            pass

    weights = root_to_leaf(raw_weights, taxa_tree)
    try:
        max_weight = max(weights.values())
    except ValueError:
        return NO_TAXA

    top_taxa = [taxa for taxa, weight in weights.items() if weight == max_weight]
    while len(top_taxa) > 1:
        lca = lowest_common_ancestor(top_taxa[0], top_taxa[1], taxa_tree)
        top_taxa = top_taxa[2:] + [lca]
    return top_taxa[0]
