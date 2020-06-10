import skbio
import pandas as pd
import pyrosetta

def fasta_to_series(path):
    '''read fasta file, return as pandas series'''
    with open(path, 'r') as f:
        data = f.read()
    data = [i.split('\n') for i in data.split('>')]
    index = [i[0] for i in data][1:] # first item is ''
    values = [''.join(i[1:]) for i in data][1:] # first item is []
    series = pd.Series(values, index=index)
    return series

def parse_mutation(m):
    '''
    m - str format "a82f" or "82f"
    not case sensitive
    returns: (num, aa)# (82, "F")
    '''
    if m != None:
        number = int(''.join([i for i in m if i.isdigit()]))
        letter = [i for i in m if i.isalpha()][-1] # in case they give 2 letters
        return number, letter.upper()
    else:
        return None

def parse_mutations(mm):
    '''
    mm: string containing more than 1 mutation (fmt:  "a82f" or "82f"), esparated by space
    returns: list of parsed mutations # [(num, aa), ...]
    '''
    return [parse_mutation(i) for i in mm.split()]

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
    for idx, (i,j) in enumerate(zip(aln_s1 ,aln_s2)):
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
    d = {i:j for i,j in enumerate(aln_s2)}
    renum = {}
    count =1
    for i in d:
        if d[i] != '-':
            renum[i] = count
            count += 1
        else:
            renum[i] = '-' # for clarity
    return renum

def map_mutations(s1,s2, mutations):
    '''
    mutations: [(pos,aa), ..] in s1
    returns mutations mapped to s2
    '''
    renum = map_sequences(s1,s2) #idx_s1:idx_s2
    print(renum)
    return [(renum[pos], aa) for pos, aa in mutations]

def mutate_pose(pose, mm):
    '''
    pose: pyrosetta pose
    mm: str of mutations like "(a)82F (f)87V" - initial letter ignored
    returns: none - inplace operation on pose
    1. make pose
    2. parse and map mutations from normal numbering to pose sequence
    3. make each mutation individually, pack side chains within 5A
    '''
    wt = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
    ACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQK\
    WERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAY\
    DENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHE\
    TTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAK\
    EDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQF\
    ALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKAEN'''.replace(' ','')

    mapped_mutations = map_mutations(wt,pose.sequence(), parse_mutations(mm))
    for pos, aa in mapped_mutations:
        print(pos, aa)
        pyrosetta.toolbox.mutate_residue(pose, pos, aa, pack_radius = 5.0)

def _test():
    series = fasta_to_series('../data/sequences/Sequences.fasta')

    wt = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
    ACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQK\
    WERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAY\
    DENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHE\
    TTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAK\
    EDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQF\
    ALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKAEN'''.replace(' ','')

    s = series.sample().item()
    d = diff(wt,s)
    print(d)
    mm = "a82f f87v l181k i263p"
    mutations = parse_mutations(mm)
    print(mutations)
    print(map_mutations(wt, s, mutations))



if __name__ == '__main__':
    _test()
