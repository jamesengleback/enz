import pyrosetta
import tools
import os

class protein():
    def __init__(self, pdb_path, seq = None):
        self.pdb_path = pdb_path
        self.clean_pdb()
        self.pdb_seq = tools.pdb_to_seq('__protein-cache__/clean.pdb')
        self.seq = seq if seq != None else tools.pdb_to_seq(self.pdb_path)
        self.map = tools.map_sequences(self.seq, self.pdb_seq)

    def mutate(self, pos, aa):
        # only changes internal sequence
        # input: self.seq position ; amino acid letter
        seq = list(self.seq) # lists are mutable
        seq[pos] = aa
        seq = ''.join(seq)
        self.seq = seq

    def mutate_pose(self, pose):
        # mutate at diff between self.seq and self.pdb_seq
        diff = tools.diff(self.pdb_seq,self.seq)
        for i in diff:
            a,b, idx = diff[i]['from'], diff[i]['to'], self.map[i]
            pyrosetta.toolbox.mutate_residue(pose, idx, b, pack_radius = 5.0)

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
        self.pose = pose

    def dump(self, path):
        if hasattr(self, 'pose'):
            self.pose.dump_pdb(path)
        else:
            self.fold()
            self.pose.dump_pdb(path)
        self.check_dump(path)

    def check_dump(self,path):
        # check if file ends with END,
        # add END if not
        # otherwise, oddt has trouble reading file
        with open(path,'r') as f:
            file = f.readlines()
        if 'END' not in file[-1]:
            # todo: check other lines for 'END'
            # blank lines not allowed
            file[-1] = file[-1].replace('\n','')
            file.append('END')
            # re-write
            with open(path,'w') as f:
                f.writelines(file)


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
        prot.mutate(i, 'A')
    prot.fold()
    print(prot.seq)

if __name__ == '__main__':
    test()
