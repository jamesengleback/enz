import pandas as pd
from tqdm import tqdm

from enz import tools
from enz import enz

bm3_wt = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRF\
DKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDM\
TRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADR\
KASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDP\
VPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERF\
ENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGG\
IPSPSTEQSAKKVRKKGC*'

def _test_Protein():
    path = 'data/clean/1jme_clean.pdb'
    wt = tools.fasta_to_series('bm3-wt.fasta')[0]
    bm3 = enz.Protein(pdb_path = 'data/clean/1jme_clean.pdb', seq = wt)
    for i in range(80,90):
        bm3.mutate_seq(i, 'A')
    bm3.refold()

def _test_vina_1():
    '''
    Test vina object with clean pdb file
    '''
    df = pd.read_csv('data/cpds/HRAC_Herbicides.csv')
    df = df.sample(3)
    vina = enz.Vina('data/clean/1jme_clean.pdb')
    for s, n in tqdm(zip(df['SMILES'],df['Name']), total = len(df)):
        scores = vina.dock(s,n)
        print(scores)

def _test_vina_2():
    '''
    test vina obj with enz.Protein
    '''
    df = pd.read_csv('data/cpds/HRAC_Herbicides.csv')
    df = df.sample(3)

    path = 'data/clean/1jme_clean.pdb'
        bm3 = enz.Protein(pdb_path = 'data/clean/3ben_clean.pdb', seq = bm3_wt)

    vina = enz.Vina(bm3)
    for s, n in tqdm(zip(df['SMILES'],df['Name']), total = len(df)):
        scores, results = vina.dock(s,n)
        print(scores)

def main():
    _test_vina_2()


if __name__ == '__main__':
    main()
