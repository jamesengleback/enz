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


def renumber(s1,s2, return_both = False):
    alns = pairwise2.align.globalxx(s1,s2,penalize_end_gaps = False, one_alignment_only=True)
    s1, s2, _,_,_ = alns[0]
    s1 = re.sub('\w-\w','',s1)
    s2 = re.sub('\w-\w','',s2)
    if return_both:
        return dict(enumerate(s1,0)), dict(enumerate(s2,0))
    else:
        return dict(enumerate(s2,0))

def test1():
    ref_seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')

    sequences = FastaToSeries('../data/sequences/Sequences.fasta')
    test_seq = sequences.sample()
    d1, d2 = renumber(ref_seq, test_seq.item(), return_both=True)

    print(test_seq.index.item())
    print(' \t','Ref ',' Match')
    for i,j in zip(d1,d2):
        print(i,'\t',d1[i],'\t',d2[j])


def test2():
    ref_seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')

    sequences = FastaToSeries('../data/sequences/Sequences.fasta')
    test_seq = sequences.sample()
    d1 = renumber(ref_seq, test_seq.item())
    print(test_seq.index.item())
    for i in d1:
        print(i,'\t',d1[i])

if __name__ == '__main__':
    test1()
