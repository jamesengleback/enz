import pandas as pd
from Bio import pairwise2
import re

def FastaToSeries(path):
    # opens fasta files, outputs series
    # full fasta string in one cell
    with open(path,'r') as file:
        data = file.read()
    data = [i.split('\n') for i in data.split('>')]
    index = [i[0] for i in data][1:] # first item is ''
    values = [''.join(i[1:]) for i in data][1:] # first item is []
    df = pd.Series(values, index=index)
    return df

def Renumber(s1,s2, return_both = False):
    # can't handle: 2 mutations in a row
    # edge cases: 1smi, non bm3 heme  sequences
    alns = pairwise2.align.globalxx(s1,s2,penalize_end_gaps = False, one_alignment_only=True)
    s1, s2, _,_,_ = alns[0]
    s1 = re.sub('\w-\w','',s1)
    s2 = re.sub('\w-\w','',s2)
    if return_both:
        # for testing
        return dict(enumerate(s1,0)), dict(enumerate(s2,0))
    else:
        return dict(enumerate(s2,0))


def test_renumber():
    ref_seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')

    sequences = FastaToSeries('../data/sequences/Sequences.fasta')

    for i, name in zip(sequences, sequences.index):
        print(name)
        d1,d2 = Renumber(ref_seq,i, return_both=True)
        for j,k in zip(d1,d2):
            if d1[j] != d2[k] and d2[k] != '-':
                print(j,d1[j], k,d2[k])

if __name__ == '__main__':
    test_renumber()
