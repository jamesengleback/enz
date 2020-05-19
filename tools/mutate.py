import argparse
import pandas as pd
from Bio import pairwise2


def FastaToDataFrame(path):
    # opens fasta files, outputs dataframe
    with open(path,'r') as file:
        data = file.read()
    data = [i.split('\n') for i in data.split('>')]
    index = [i[0] for i in data][1:] # first item is ''
    values = [list(''.join(i[1:])) for i in data][1:] # first item is []
    df = pd.DataFrame(values, index=index)
    return df

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

def AlignAgainstReferenceSequence(query_seq):
    reference_seq='''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQ\
    ALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRL\
    TLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADR\
    KASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVL\
    VDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEE\
    FRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKA\
    KSKKIPLGGIPSPSTEQSAKKVRKKGC*'''
    alignments = pairwise2.align.globalxx(reference_seq, query_seq) # makes several alignments
    all_scores = [i[2] for i in alignments]
    top_scores_indexes = [i for i,j in enumerate(all_scores) if j ==max(all_scores)]# list of indexes where score is max()
    top_alignment = alignments[top_scores_indexes[0]] # Highest scoring first

    dictionary = {'aln_reference_seq': top_alignment[0],
    'aln_query_seq':top_alignment[1],
        'aln_score': top_alignment[2]}
    return dictionary

def ResidueConservation(seq1, seq2):
    ## takes two aligned sequences
    # probably won't work if  they're different lengths
    alignment_df = pd.DataFrame([list(seq1), list(seq2)])
    alignment_df.replace(' ','-', inplace = True)
    conserved_residue_count = 0
    for i in alignment_df:
        if len(alignment_df[i].unique()) <2:
            conserved_residue_count += 1
    frac_conserved = conserved_residue_count/len(seq1)
    return frac_conserved



def FindMutations(query_seq):
    reference_seq='''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQ\
    ALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRL'''

    alignments = pairwise2.align.localxx(reference_seq, query_seq) # makes several alignments
    all_scores = [i[2] for i in alignments]
    top_scores_indexes = [i for i,j in enumerate(all_scores) if j ==max(all_scores)]# list of indexes where score is max()
    top_alignment = alignments[top_scores_indexes[0]] # Highest scoring first

    dictionary = {'aln_reference_seq': top_alignment[0].replace(' ','-'),
    'aln_query_seq':top_alignment[1].replace(' ','-'),
        'aln_score': top_alignment[2]}

    mutations = [{'ref':j[0],'pos':i,'query':j[1]} for i, j in enumerate(zip(dictionary['aln_reference_seq'],\
    dictionary['aln_query_seq'])) if j[0] != j[1] ] # list of dictionarys
    # trim down
    mutations = [i for i in mutations if i['pos'] > 20 and i['pos'] < 400]
    print(mutations)


def main(args):
    BM3s = FastaToSeries(args.sequences)
    mutant_id = [i for i in BM3s.index if '1p0x' in i][0] # should be one item
    for i in BM3s:
        mutant_sequence = i

        alignment = AlignAgainstReferenceSequence(mutant_sequence)
        frac_conserved = ResidueConservation(alignment['aln_reference_seq'], alignment['aln_query_seq'])
        print(f'Percentage Conservation: {round(100 * frac_conserved,2)} %')

        FindMutations(mutant_sequence)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--sequences', help = 'MSA file')
    args = parser.parse_args()
    main(args)
