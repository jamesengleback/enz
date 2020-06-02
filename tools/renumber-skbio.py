import pandas as pd
import skbio

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

def Renumber(query_seq, return_both=False):
    # canonical numbering
    ref_seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')
    # sequences need to be in the skbio.Protein object
    ref_seq = skbio.Protein(ref_seq)
    query_seq = skbio.Protein(query_seq) # item unpacks the string
    # align - slow, gives efficiency warning - skbio devs working on this
    # check for updates in future
    print('\x1b[31m') #red
    aln = skbio.alignment.global_pairwise_align_protein(ref_seq,query_seq)
    print('\x1b[0m') # reset colors
    # unpack aligned sequences
    aln_ref = ''.join([i.decode() for i in aln[0].loc[0].values])
    aln_query = ''.join([i.decode() for i in aln[0].loc[1].values])
    if return_both:
        return dict(enumerate(aln_query,0)), dict(enumerate(aln_ref,0))
    else:
        return dict(enumerate(aln_query,0))

def IndexMutations(query_seq):
    # returns dictionary of renumbered mutations
    # get renumbered dictionaries of query and mutant
    d_query, d_ref = Renumber(query_seq, return_both=True)
    # find differences (not including gaps)
    mutations = {}
    for i, j in zip(d_query, d_ref):
        query_aa = d_query[i]
        ref_aa = d_ref[j]
        if query_aa != ref_aa and query_aa != '-' and ref_aa != '-':
            mutations[i] = {'from':ref_aa, 'to':query_aa}
    return mutations

def MapMutations(query_seq):
    # 1. renumber query_seq - map to original
    # 2. index mutations - use same alignment
    # 3. return mutation positions in
    d_query, d_ref = Renumber(query_seq, return_both=True)
    tmp = {i:d_query[i] for i in d_query if d_query[i] != '-'}
    mapping = {j:i for i,j in enumerate(tmp, 1)} # original:renumbered, start count at 1
    mutations = {}
    for i, j in zip(d_query, d_ref):
        query_aa = d_query[i]
        ref_aa = d_ref[j]
        if query_aa != ref_aa and query_aa != '-' and ref_aa != '-':
            mutations[i] = {'from':ref_aa, 'to':query_aa}

    mapped_mutations = {mapping[i]:mutations[i]['to'] for i in mutations}
    return mapped_mutations


def test():
    sequences = FastaToSeries('../data/sequences/Sequences.fasta')
    query_seq = sequences.sample().item()
    #renumbered_seq = Renumber(query_seq)
    #print(renumbered_seq)
    #mutations = IndexMutations(query_seq)
    #print(mutations)
    mapping = MapMutations(query_seq)
    print(mapping)


if __name__ == '__main__':
    test()
