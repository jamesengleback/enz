#import sys
#sys.path.append('../enz')
import enz

seq = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
ACDESRFDKNLSQAWKFVRDFAGDGLVTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKSQRANPDDP\
AYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIA\
GHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFS\
LYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRAC\
IGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK'

def test_protein():
    pdb = '../data/1jme.pdb'
    p = enz.protein(pdb, seq = seq, key_sites = [82,87,188,330])
    print(p.df)
    p.mutate(24,'A')
    p.refold()
    print(p.KEY_SITES_DICT)
    print(p.seq)
    print(p.df)
    print(p.docking_results)
    print(p.dock('CCCCCC=O'))
    p.save('tmp.pdb')
    p.save_docking_results('results')


def test_vina_pdb():
    pdb = '../data/1jme.pdb'
    v = enz.vina(pdb)
    v.dock('CCCC=O')
    print(v.df)
    v.save('results2')

def test_vina_enz_protein():
    pdb = '../data/4key.pdb'
    p = enz.protein(pdb, seq = seq, key_sites = [82,87,188,330,75,181,49,51], cofactors = ['HEM'])
    print(p.dock('CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3'))

def pdb_cleaning():
    pdb = '../data/4key.pdb'
    p = enz.protein(pdb, seq = seq, key_sites = [82,87,188,330,75,181,49,51], cofactors = ['HEM'])
    print(p.df)

def main():
    pdb_cleaning()



if __name__ == '__main__':
    main()
