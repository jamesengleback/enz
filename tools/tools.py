import skbio
import pandas as pd

def fasta_to_series(path):
    with open(path, 'r') as f:
        data = f.read()
    data = [i.split('\n') for i in data.split('>')]
    index = [i[0] for i in data][1:] # first item is ''
    values = [''.join(i[1:]) for i in data][1:] # first item is []
    series = pd.Series(values, index=index)
    return series

def aln(s1,s2):
    '''align s1 and s2, extract and return aligned sequences'''
    aln = skbio.alignment.global_pairwise_align_protein(skbio.Protein(s1),skbio.Protein(s2))[0]
    aln_s1 = ''.join([i.decode() for i in aln[0].values])
    aln_s2 = ''.join([i.decode() for i in aln.loc[1].values])
    return aln_s1, aln_s2

def diff(s1,s2):
    '''
    return mutations in s2, relative to s1;
    ignores gaps (indels);
    returns dict'''
    aln_s1 ,aln_s2 = aln(s1,s2)
    m = {}
    for idx, (i,j) in enumerate(zip(aln_s1 ,aln_s2),1):
        if i != j and i != '-' and j != '-':
            m[idx] = {'from':i,'to':j}
    return m


def mutate_sequence(seq, pos,aa):
    # single mutation
    # unnecessary
    s = dict(enumerate(seq,1))
    s[pos] = aa
    return ''.join(s.values())

def map_sequences(s1,s2):
    '''
    s1 & s2 - ref & query sequences - str
    returns {idx:idx} mapping for s1:s2
    '''
    aln_s1 ,aln_s2= aln(s1,s2)
    d = {i:j for i,j in enumerate(aln_s2,1)}
    renum = {}
    count =1
    for i in d:
        if d[i] != '-':
            renum[i] = count
            count += 1
    return renum

def _test():
    s1 = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL'
    s2 = 'MTIPEMPQPKPVQALMKGEIFKFEPPGRPTYLSSQRL'
    x = map_sequences(s1,s2)
    print(x)
    series = fasta_to_series('../data/sequences/Sequences.fasta')

    wt = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
    ACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQK\
    WERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAY\
    DENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHE\
    TTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAK\
    EDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQF\
    ALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKAEN'''.replace(' ','')

    from tqdm import tqdm
    idx =[]
    mutations = []
    for i in tqdm(series.index):
        d = diff(wt,series[i])
        fmt_mutations = [f"{d[j]['from']}{j}{d[j]['to']}" for j in d]
        pdb_id = i.split('_')[0]
        mutations.append(fmt_mutations)
        idx.append(pdb_id)
    pd.Series(mutations, index = idx).to_csv()


if __name__ == '__main__':
    _test()
