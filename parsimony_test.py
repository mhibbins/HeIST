
from Bio import Phylo
from io import StringIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import * 

test_newick_phy = Phylo.read(StringIO('(1:6,(2:3,((6:1.25,5:1.25):0.25,(4:1,3:1):0.5):1.5):3)'), 'newick')

a = SeqRecord(Seq('A', generic_dna), id='1')
b = SeqRecord(Seq('T', generic_dna), id='2')
c = SeqRecord(Seq('T', generic_dna), id='3')
d = SeqRecord(Seq('T', generic_dna), id='4')
e = SeqRecord(Seq('A', generic_dna), id='5')
f = SeqRecord(Seq('T', generic_dna), id='6')

test_pattern = MultipleSeqAlignment([a,b,c,d,e,f])


def get_min_mutations(tree, traitpattern):
    '''
    Gets the minimum number of mutations required
    to explain the trait pattern without hemiplasy;
    ie. the parsimony score. Takes a newick tree 
    and trait pattern in alignment form
    '''
    scorer = ParsimonyScorer()
    
    return scorer.get_score(tree, traitpattern)
    #return ParsimonyScorer.get_score(tree, traitpattern)    
    

test = get_min_mutations(test_newick_phy, test_pattern)
print(test)
