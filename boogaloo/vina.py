#!/home/james/miniconda3/envs/enz2/bin/python
import sys
import os
import subprocess
import tempfile
from openbabel import pybel

'''
todo: 
    residue select - for box - requrires a pdb parser
    save output
    score - parse vina output:

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         11.8      0.000      0.000
   2         12.8      1.122      1.639
   3         13.8      2.289      4.663
   4         13.9      1.860      2.665
   5         14.1      2.627      4.719
   6         14.6      1.578      2.352

'''

class Vina:
    '''
    runs vina as subprocess
    accepts pdb and smiles as input
    todo: residue numbering and center finding
    '''
    def __init__(self, receptor, ligand):
        self.EX = '/home/james/src/autodock_vina_1_1_2_linux_x86/bin/vina'
        self.CACHE = tempfile.mkdtemp(prefix='boogaloo')
        self.receptor = self.prep_receptor(receptor)
        self.ligand = self.prep_ligand(ligand)

        self.args = {'--receptor':self.receptor,
                    '--ligand':self.ligand,
                    '--center_x':10,
                    '--center_y':10,
                    '--center_z':10,
                    '--size_x':10,
                    '--size_y':10,
                    '--size_z':10,
                    '--out':None,
                    '--exhaustiveness':None}

    def prep_receptor(self, receptor):
        m = list(pybel.readfile('pdb',receptor)) 
        # pybel.readfile returns generator of mols in file
        # hopefully just 1
        assert len(m) == 1
        m = m[0]
        m.addh()
        savepath = os.path.join(self.CACHE, 'receptor.pbqt')
        m.write('pdbqt',savepath, opt={'r':True}, overwrite=True)
        # opt r save pdbqt as rigid, prevents parse errors
        return savepath

    def prep_ligand(self, smiles):
        m = pybel.readstring('smi',smiles)
        m.addh() # redundant - included in make3D()
        m.make3D()
        savepath = os.path.join(self.CACHE, 'ligand.pdbqt')
        m.write('pdbqt', savepath)
        return savepath

    def dock(self):
        args = [self.EX]
        for i in self.args:
            if self.args[i] != None:
                args.append(i)
                args.append(str(self.args[i]))
        #### issue - pasring receptor - ROOT - inapproproate tag
        ### appantly obabel -r sets receptor as rigid and fizes this
        print(args)
        print(subprocess.Popen(args))


EX = '/home/james/src/autodock_vina_1_1_2_linux_x86/bin/vina'
prot = 'test_data/4key.pdb'
def main():

    vina = Vina(prot, 'CCCCCCC=O')
    vina.dock()
    #subprocess.check_output([EX, '--receptor', prot, '--ligand',lig])

if __name__ == '__main__':
    main()
