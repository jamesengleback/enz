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

class mol:
    def __init__(self, path):
        pass
    @property
    def df(self):
        pass
    def save(self, path):
        pass

class protein:
    # front end
    def __init__(self, pdb_path, seq = None, cofactors = [], key_sites = []):
        self.pdb_path = pdb_path
        self.CACHE = tempfile.mkdtemp()
        self.key_sites = key_sites
        self.cofactors = cofactors
        self.struc = pdb_fns.clean_pdb(pdb_path = self.pdb_path,
                    save_path = os.path.join(self.CACHE, 'clean.pdb'),
                    cofactors = self.cofactors)
        self.pdb_seq = pdb_fns.get_seq(self.struc)
        self.seq = self.pdb_seq if seq == None else seq
    @property
    def df(self):
        pass
    def mutate(self, position, aa):
        seq = list(self.seq)
        seq[position] = aa
        self.seq = ''.join(seq)
    def refold(self):
        aln1, aln2 = aln(self.pdb_seq, self.seq)
        mutations = diff(aln1, aln2)
        global PYROSETTA_INIT
        if  not PYROSETTA_INIT:
            pyrosetta_init(silent=True)
            PYROSETTA_INIT = True
        self.pose = pose_from_pdb(self.struc)
        for i in mutations:
            mutate_residue(self.pose, i, mutations[i]['to'], pack_radius = 5.0)
        self.pose.dump_pdb(self.struc)
    def save(self, path):
        # copy self.stuc -> path
        shutil.copyfile(self.struc, path)
    def dock(self, smiles, save_path = None):
        df = vina.dock(self.struc,
                    smiles,
                    save_path = save_path,
                    cofactors = self.cofactors,
                    target_residues = self.key_sites)
        return df

class pdb_fns:
    def clean_pdb(pdb_path, save_path, cofactors = [], chain_selection = 'A'):
        structure = PandasPdb().read_pdb(pdb_path)
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
    def get_seq(pdb_path):
        structure = PandasPdb().read_pdb(pdb_path)
        sequences = structure.amino3to1() # cols = ['chain_id', 'residue_name']
        seqs = [''.join(sequences.loc[sequences['chain_id'] == i,'residue_name'].to_list()) for i in sequences['chain_id'].unique()]
        return seqs[0] if len(seqs) == 1 else seqs
    def draw_box(pdb_path, key_sites):
        receptor = PandasPdb().read_pdb(pdb_path)
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

class vina:
    def dock(receptor_pdb,
            smiles,
            save_path = None,
            cofactors = [],
            target_residues = [],
            exhaustiveness=8,
            vina_executable = find_executable('vina'),
            vina_split_executable = find_executable('vina_split')):
        CACHE = tempfile.mkdtemp()
        raw_vina_results = os.path.join(CACHE, 'vina.result')
        # todo : if not clean
        clean_receptor_pdb = pdb_fns.clean_pdb(receptor_pdb, os.path.join(CACHE, f'{os.path.basename(receptor_pdb)}.clean'), cofactors = cofactors)
        receptor_pdbqt = obabel_fns.pdb_to_pdbqt(clean_receptor_pdb, os.path.join(CACHE,'receptor.pdbqt'))
        ligand_pdbqt = obabel_fns.smiles_to_pdbqt(smiles, os.path.join(CACHE,'ligand.pdbqt'))
        args = {'--receptor':receptor_pdbqt,
                    '--ligand':ligand_pdbqt,
                    '--out':raw_vina_results,
                    '--exhaustiveness':exhaustiveness}
        args.update(pdb_fns.draw_box(clean_receptor_pdb, target_residues))
        args_list_vina = [vina_executable]
        for i in args:
            args_list_vina.append(i)
            args_list_vina.append(str(args[i]))
        # execute
        p1 = subprocess.check_output(args_list_vina)
        docking_scores = vina.extract_scores(p1.decode())
        if save_path != None:
            vina.split(raw_vina_results,
                        save_path,
                        vina_split_executable)
        return docking_scores

    def split(vina_output, save_path, vina_split_executable):
        # vina_split
        args_list_vina_split = [vina_split_executable, '--input', outpath]
        p2 = subprocess.Popen(args_list_vina_split, stdout=subprocess.DEVNULL)
        p2.wait()
        # convert pdbqt files
        docking_poses_pdbqt = [os.path.join(CACHE,i) for i in os.listdir(CACHE) if 'vina_ligand_' in i]
        docking_poses_pdb = [save_pose(i, os.path.join(CACHE, f'pose_{j}.pdb')) for j,i in enumerate(docking_poses_pdbqt)]
        # move to savepath
        os.makedirs(save_path, exist_ok=True)
        # remove existing contents
        for i in os.listdir(save_path):
            os.remove(os.path.join(save_path,i))
        for i in docking_poses_pdb:
            shutil.copyfile(i, os.path.join(save_path, os.path.basename(i))) # copy to inside savepath dir
        shutil.copyfile(clean_receptor_pdb, os.path.join(save_path, os.path.basename(clean_receptor_pdb)))
        docking_scores.to_csv(os.path.join(save_path, 'docking-scores.csv')) # and save scores

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

def aln(s1,s2):
    return nw.global_align(s1,s2)

def diff(s1,s2):
    return {i:{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s1,s2)) if x != y and x != '-' and y != '-'}

def test():
    path = '../test/4KEY.pdb'
    p = protein(path, cofactors = 'HEM', key_sites = [87,82,400,330,263,188,49,51])
    p.mutate(70,'A')
    p.refold()
    p.save('test.pdb')
    print(p.dock('CCCCCCCCCCCC=O', 'results'))

if __name__ == '__main__':
    test()
