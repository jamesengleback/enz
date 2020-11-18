import tempfile
import os
import shutil
from distutils.spawn import find_executable
import subprocess
import re

import pandas as pd

from biopandas.pdb import PandasPdb
import nwalign3 as nw
from openbabel import pybel

from pyrosetta import pose_from_pdb
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue

PYROSETTA_INIT = False
pybel.ob.obErrorLog.SetOutputLevel(0)

'''
todo 
  - results object could be simpler
'''
class mol:
    '''
    protein & results poses inherit from this class
    '''
    def __init__(self, struc):
        self.struc = struc
    @property
    def df(self):
        data = PandasPdb().read_pdb(self.struc)
        return data.df['ATOM'].append(data.df['HETATM']) # all atoms
    def save(self, save_path):
        shutil.copyfile(self.struc, save_path)


class protein(mol):
    '''
    stores clean protein sequence & structure
    cleaning removes waters and hetatm residues not in cofactors
    used for mutant structre prediction and small molecule docking
    example:
    >>> import enz
    >>> wt = 'MTIKEM...'
    >>> 
    >>> for i in [75,87,330, 263]:
    >>>    p = enz.protein('enz/data/4key.pdb', seq=wt)
    >>>    p.mutate(i, 'A') # alanine scan
    >>>    r = p.dock('CCCCCCC=O', target_residues=[82,87,330,400,51]) # specify docking site
    >>>    r.save(f'bm3_{p.seq[i]}iA') # e.g. bm3_F87A

    '''
    def __init__(self, struc, seq = None, cofactors = [], key_sites = []):
        super().__init__(struc)
        self.key_sites = key_sites
        self.cofactors = cofactors
        self.CACHE = tempfile.mkdtemp()
        self.struc = pdb_fns.clean_pdb(struc = self.struc,
                    save_path = os.path.join(self.CACHE, 'clean.pdb'),
                    cofactors = self.cofactors)
        self.pdb_seq = pdb_fns.get_seq(self.struc)
        self.seq = self.pdb_seq if seq == None else seq

    def mutate(self, position, aa):
        seq = list(self.seq)
        seq[position] = aa
        self.seq = ''.join(seq)
    
    def refold(self):
        aln1, aln2 = utils.aln(self.pdb_seq, self.seq)
        mutations = utils.diff(aln1, aln2)
        global PYROSETTA_INIT
        if  not PYROSETTA_INIT:
            pyrosetta_init(silent=True)
            PYROSETTA_INIT = True
        self.pose = pose_from_pdb(self.struc) # pyrosetta
        for i in mutations:
            mutate_residue(self.pose, i, mutations[i]['to'], pack_radius = 5.0)
        self.pose.dump_pdb(self.struc)
    
    def dock(self, smiles, save_path = None, target_residues = None):
        if target_residues == None:
            target_residues = self.key_sites
        results = vina.dock(self.struc,
                    smiles,
                    save_path = save_path,
                    cofactors = self.cofactors,
                    target_residues = target_residues)
        results.receptor = self
        return results

class pdb_fns:
    def clean_pdb(struc, save_path, cofactors = [], chain_selection = 'A'):
        structure = PandasPdb().read_pdb(struc)
        atoms = structure.df['ATOM'].copy()
        hetatms = structure.df['HETATM'].copy()
        atoms = atoms.loc[atoms['chain_id'] == chain_selection,:]
        hetatms = hetatms.loc[hetatms['chain_id'] == chain_selection,:]
        het_garbage = [i for i in hetatms['residue_name'].unique() if i not in cofactors]
        hetatms = hetatms.loc[hetatms['residue_name'].isin(het_garbage) == False,:]
        structure.df['ATOM'] = atoms
        structure.df['HETATM'] = hetatms
        structure.to_pdb(save_path)
        return save_path

    def get_seq(struc):
        structure = PandasPdb().read_pdb(struc)
        sequences = structure.amino3to1() # cols = ['chain_id', 'residue_name']
        seqs = [''.join(sequences.loc[sequences['chain_id'] == i,'residue_name'].to_list()) for i in sequences['chain_id'].unique()]
        return seqs[0] if len(seqs) == 1 else seqs

    def draw_box(struc, key_sites):
        receptor = PandasPdb().read_pdb(struc)
        df = receptor.df['ATOM']
        target_site = df.loc[df['residue_number'].isin(key_sites),:]
        coords = target_site.loc[:,['x_coord','y_coord','z_coord']]
        center = coords.mean(axis=0)
        sizes = (coords.max(axis=0) - coords.min(axis=0)) * 1.2
        box = {'--center_x':center['x_coord'],
                '--center_y':center['y_coord'],
                '--center_z':center['z_coord'],
                '--size_x':sizes['x_coord'],
                '--size_y':sizes['y_coord'],
                '--size_z':sizes['z_coord']}
        return box

class obabel_fns:
    def pdb_to_pdbqt(pdb, save_path):
        m = list(pybel.readfile('pdb',pdb))
        assert len(m) == 1
        m = m[0]
        m.addh()
        m.write('pdbqt', save_path, opt={'r':True}, overwrite=True) # opt:r = rigid - less errors?? - revisit this
        return save_path

    def smiles_to_pdbqt(smiles, save_path):
        m = pybel.readstring('smi',smiles)
        m.addh()
        m.make3D()
        m.write('pdbqt',save_path)
        return save_path

    def pdbqt_to_pdb(pdbqt, save_path):
        m = list(pybel.readfile('pdbqt',pdbqt))
        assert len(m) == 1
        m = m[0]
        # already 3d
        m.write('pdb', save_path)


class vina:
    def dock(receptor_pdb,
            smiles,
            save_path = None,
            cofactors = [],
            target_residues = [],
            exhaustiveness=8,
            vina_executable = find_executable('vina'),
            vina_split_executable = find_executable('vina_split')):
        # check there's a box
        if target_residues == []:
            raise Exception('no target residues selected')

        CACHE = tempfile.mkdtemp()
        raw_vina_results = os.path.join(CACHE, 'vina.result')
        # todo : if not clean  - check if dock from protein object
        clean_receptor_pdb = pdb_fns.clean_pdb(receptor_pdb, os.path.join(CACHE, f'{os.path.basename(receptor_pdb)}.clean'),
                                                                            cofactors = cofactors)
        receptor_pdbqt = obabel_fns.pdb_to_pdbqt(clean_receptor_pdb, os.path.join(CACHE,'receptor.pdbqt'))
        ligand_pdbqt = obabel_fns.smiles_to_pdbqt(smiles, os.path.join(CACHE,'ligand.pdbqt'))
        args = {'--receptor':receptor_pdbqt,
                    '--ligand':ligand_pdbqt,
                    '--out':raw_vina_results,
                    '--exhaustiveness':exhaustiveness}
        args.update(pdb_fns.draw_box(clean_receptor_pdb, target_residues)) # add box dims to args

        args_list_vina = [vina_executable]
        for i in args:
            args_list_vina.append(i)
            args_list_vina.append(str(args[i]))

        # execute
        p1 = subprocess.check_output(args_list_vina)
        
        # create results object
        docking_scores = vina.extract_scores(p1.decode())
        poses = vina.vina_split(raw_vina_results, vina_split_executable)
        results = vina.results([os.path.join(poses, i) for i in os.listdir(poses)], docking_scores)
        
        return results

    def vina_split(raw_vina_results, vina_split_executable):
        # vina_split
        args_list_vina_split = [vina_split_executable, '--input', raw_vina_results]
        p = subprocess.Popen(args_list_vina_split, stdout=subprocess.DEVNULL)
        p.wait() # outputs odbqt files
        results_dir = os.path.dirname(raw_vina_results)
        poses = [os.path.join(results_dir, i) for i in os.listdir(results_dir) if 'vina.result_ligand' in i]
        clean_results = os.path.join(results_dir, 'pose_pdbs')
        os.mkdir(clean_results)
        for i in poses:
            save_path = os.path.basename(i).replace('pdbqt','pdb')
            obabel_fns.pdbqt_to_pdb(i, os.path.join(clean_results, save_path))
        return clean_results # path to foder containing pdb of poses

    def extract_scores(text):
        # extract scores from vina output
        text = text.split('\n')
        table_start = ['---+--' in i for i in text].index(True) + 1
        is_all_ints = lambda l : sum([re.search('-?\d+', i) is not None for i in l])
        table = []
        for row in text[table_start:]:
            items = row.split()
            if len(items) == 4 and is_all_ints(items):
                table.append(dict(zip(['mode','affinity (kcal/mol)', 'dist from best mode - rmsd - ub','dist from best mode - lb'], items)))
        return pd.DataFrame(table)

    class results:
        '''
        poses & score df
        '''
        def __init__(self, poses, scores):
            self.poses = {int(re.findall('\d+',os.path.basename(i))[0])\
                    :mol(i) for i in poses}
            self.scores = scores.astype(float)
            self.dictionary = {os.path.basename(i.struc):{'mol':i, 'affinity':j} for i,j in zip(self.poses.values(), self.scores['affinity (kcal/mol)'])}
        def save(self, save_path):
            os.makedirs(save_path, exist_ok = True)
            self.scores.to_csv(os.path.join(save_path, 'scores.csv'))
            for i in self.poses:
                pose_i = self.poses[i]
                pose_i.save(os.path.join(save_path, os.path.basename(pose_i.struc)))
            self.receptor.save(os.path.join(save_path, 'receptor.pdb'))

class utils:
    def aln(s1, s2):
        aln1, aln2 = nw.global_align(s1,s2)
        return aln1, aln2

    def diff(s1,s2):
        return {i:{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s1,s2)) if x != y and x != '-' and y != '-'}


def test():
    bmw_wt = 'TIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK*'
    path = '../data/4key.pdb'
    smiles = 'CCCCCCCC=O'
    p = protein(path, cofactors = ['HEM'], seq=bmw_wt)
    #p.mutate(82,'F')
    #p.refold()
    r = p.dock(smiles, target_residues = [82,87,400,188,181,263])
    print(r)
    r.save('test')
if __name__ == '__main__':
    test()

