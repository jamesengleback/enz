import pyrosetta
import tools
import os

class protein():
    def __init__(self, pdb_path, seq = None):
        self.pdb_path = pdb_path
        self.clean_pdb()
        self.pdb_seq = tools.pdb_to_seq('__protein-cache__/clean.pdb')
        self.seq = seq if seq != None else tools.pdb_to_seq(self.pdb)
        self.map = tools.map_sequences(self.seq, self.pdb_seq)

    def mutate_seq(self, pos, aa):
        # input: self.seq position ; amino acid letter
        seq = list(self.seq) # lists are mutable
        seq[pos] = aa
        seq = ''.join(seq)
        self.seq = seq

    def mutate_pose(self, pose):
        # mutate at diff between self.seq and self.pdb_seq
        for i in self.map:
            correspond = self.map[i]
            if correspond != '-':
                a = self.seq[i] # true seq
                b = self.pdb_seq[correspond] # pdb seq
                if a != b and a != '-' and b != '-':
                    print(a,b)
                    #pyrosetta.toolbox.mutate_residue(pose, correspond, a, pack_radius = 5.0)

    def clean_pdb(self):
        pdb = tools.clean_pdb(self.pdb_path)
        cache = os.makedirs('__protein-cache__', exist_ok=True)
        pdb.to_pdb(os.path.join('__protein-cache__', 'clean.pdb'))

    def fold(self):
        # n decoys?
        # n cpus?
        pyrosetta.init(silent=True)
        pose = pyrosetta.pose_from_pdb('__protein-cache__/clean.pdb')
        # mutate
        self.mutate_pose(pose)
        # loop remodel

    def dump(self):
        pass

def test():
    wt = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
    ACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQK\
    WERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAY\
    DENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHE\
    TTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAK\
    EDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQF\
    ALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKAEN'''.replace(' ','')

    prot = protein(pdb_path = '../data/clean/1jme_clean.pdb', seq = wt)
    for i in range(80,90):
        prot.mutate_seq(i, 'A')
    prot.fold()
    print(prot.seq)

if __name__ == '__main__':
    test()
