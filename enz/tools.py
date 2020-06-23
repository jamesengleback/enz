import skbio
import pandas as pd
from biopandas.pdb import PandasPdb
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

def pdb_to_seq(pdb):
    pdb = PandasPdb().read_pdb(pdb)
    df = pdb.amino3to1() # cols: 'chain_id', 'residue_name'
    assert len(df['chain_id'].unique()) == 1
    seq = ''.join([i for i in df['residue_name']])
    return seq

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
    return [(renum[pos], aa) for pos, aa in mutations]

def mutate_pose(pose, seq, mm): ###########################
    '''
    pose: pyrosetta pose
    seq: string - sequence
    mm: list of mutations - may need parsing with parse_mutation first
    1. make pose
    2. parse and map mutations from normal numbering to pose sequence
    3. make each mutation individually, pack side chains within 5A
    '''

    mapped_mutations = map_mutations(seq,pose.sequence(), mm)
    for pos, aa in mapped_mutations:
        print(pos, aa)
        pyrosetta.toolbox.mutate_residue(pose, pos, aa, pack_radius = 5.0)


def clean_pdb(pdb_path):
    data = PandasPdb().read_pdb(pdb_path)
    singleChainData = SingleOutChain(data)
    strippedDown = StripHetAtoms(singleChainData)
    return singleChainData


def SingleOutChain(pandas_pdb):
    # return complete pandas pdb object
    unique_chains = pandas_pdb.df['ATOM']['chain_id'].unique()
    chain_1 = unique_chains[0]
    pandas_pdb.df['ATOM'] = pandas_pdb.df['ATOM'].loc[pandas_pdb.df['ATOM']['chain_id'] == chain_1]
    return pandas_pdb

def StripHetAtoms(pandas_pdb):
    # strips water and ligands, except the heme
    unique_chains = pandas_pdb.df['ATOM']['chain_id'].unique()
    chain_1 = unique_chains[0]

    pandas_pdb.df['HETATM'] = pandas_pdb.df['HETATM'].loc[pandas_pdb.df['HETATM']['residue_name'] == 'HEM']
    pandas_pdb.df['HETATM'] = pandas_pdb.df['HETATM'].loc[pandas_pdb.df['HETATM']['chain_id'] == chain_1]
    return pandas_pdb

    
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
    pyrosetta.init()
    pose = pyrosetta.pose_from_pdb('../data/clean/1jme_clean.pdb')
    mutate_pose(pose, wt, mutations)
    print(pose.sequence())



if __name__ == '__main__':
    _test()
