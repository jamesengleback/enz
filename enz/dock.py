import os
import numpy as np
import pandas as pd
from tqdm import  tqdm
import multiprocessing

import oddt
from oddt import docking
from oddt.scoring.functions import NNScore

import protein

class Vina():
    def __init__(self, protein = None, screen=None):
        self.cache = '__vina-cache__'
        os.makedirs(self.cache, exist_ok=True)
        self.receptor = self.process_protein(protein)
        if screen != None:
            self.screen = self.read_screen(screen)
    def process_protein(self,prot):
        # type: pdb path, enz.protein.protein
        # return oddt object
        if isinstance(prot, protein.protein):
            oddt_protein = self.read_protein(prot)
        else:
            # assuming protein is pdb path
            oddt_protein = self.read_pdb(prot)
        return oddt_protein

    def read_protein(self,protein):
        # if enz.protein: dump in cache
        path = os.path.join(self.cache, 'receptor.pdb')
        protein.dump(path) # save as pdb
        oddt_protein = self.read_pdb(path)
        return oddt_protein

    def read_pdb(self,path):
        print(path)
        oddt_protein = next(oddt.toolkit.readfile('pdb', path))
        oddt_protein.protein = True # register that it's a protein
        oddt_protein.addh()
        return oddt_protein
    def define_active_site(self):
        pass
    def read_screen(self,screen):
        # expected type: screen == pd.DataFrame
        # columns = smiles, names
        smiles_col = [i for i in screen.columns if 'smiles' in i.lower()]
        assert len(smiles_col) == 1
        name_col = [i for i in screen.columns if 'smiles' not in i.lower()]
        assert len(name_col) == 1
        smiles = screen[smiles_col[0]]
        names = screen[name_col[0]]
        self.screen = {'smiles':smiles,'names':names}

    def read_smiles(self,smiles):
        mol = oddt.toolkit.readstring('smi',smiles)
        mol.addh()
        mol.make3D()
        return mol

    def autodock_score(self, results, vina):
        scores = pd.concat([pd.Series(i.data) for i in vina.score(results)],
        axis=1,
        join='outer').T
        return scores
    def init_vina(self, ncpu):
        if not hasattr(self,'vina'):
            self.vina = docking.autodock_vina(self.receptor, n_cpu = ncpu) # oddt.docking
        # todo: check on ncpus
    def dock(self, ncpu = multiprocessing.cpu_count() -1):
        if not hasattr(self,'vina'):
            self.init_vina(ncpu)
        # scores
        df = pd.DataFrame([],columns=['vina_affinity', 'vina_rmsd_lb', 'vina_rmsd_ub', 'vina_rmsd_input',
           'vina_rmsd_input_min', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion',
           'vina_hydrophobic', 'vina_hydrogen', 'SMILES Atom Order', 'cpd'])

        scores_path = os.path.join(self.cache,'autodock_scores.csv')

        if not os.path.exists(scores_path):
            df.to_csv(scores_path)

        for name,smiles in tqdm(zip(self.screen['names'], self.screen['smiles']),total=len(self.screen['smiles'])):
            results = self.vina.dock(self.read_smiles(smiles))
            scores = self.autodock_score(results, self.vina)
            scores['cpd'] = name
            print(scores)
            scores.fillna('-')
            scores.to_csv(scores_path, mode='a', header=None)



def read_smiles(smiles):
    mol = oddt.toolkit.readstring('smi',smiles)
    mol.addh()
    mol.make3D()
    return mol

def save_results(results, dirname=None, parent_dir = None):
    # mkdir
    if dirname==None:
        DIR = 'tmp'
    else:
        DIR=dirname
    if parent_dir != None:
        DIR = os.path.join(parent_dir,DIR)
    os.makedirs(DIR, exist_ok=True) #overwrite

    for i,mol in enumerate(results):
        name = f'{i}.pdb'
        path = os.path.join(DIR,name)
        mol.write('pdb',path,overwrite=True)

def make_protein(path):
    protein = next(oddt.toolkit.readfile('pdb', path))
    protein.protein = True # register that it's a protein
    protein.addh()
    return protein

def autodock_score(results, vina):
    scores = pd.concat([pd.Series(i.data) for i in vina.score(results)],
    axis=1,
    join='outer').T
    return scores

def dock_single(protein, oddt_mol):
    vina = docking.autodock_vina(protein)
    print('Docking..')
    results = vina.dock(oddt_mol)
    print('Done')
    save_results(results)
    SCORES = autodock_score(results, vina)
    return SCORES

def dock_multi(protein, smiles, names, target=None):
    # protein must be oddt type
    if target != None:
        outdir = 'docking_results'
    else:
        outdir = os.path.join(target,'docking_results')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    scores_path = os.path.join(outdir,'autodock_scores.csv')

    df = pd.DataFrame([],columns=['vina_affinity', 'vina_rmsd_lb', 'vina_rmsd_ub', 'vina_rmsd_input',
       'vina_rmsd_input_min', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion',
       'vina_hydrophobic', 'vina_hydrogen', 'SMILES Atom Order', 'cpd', 'tgt'])

    if not os.path.exists(scores_path):
        df.to_csv(scores_path, index=None)

    vina = docking.autodock_vina(protein)

    for i,smi in tqdm(zip(names, smiles), total=len(smiles)):
        mol = read_smiles(smi)
        results = vina.dock(mol)
        save_results(results, i, parent_dir=outdir)
        scores = autodock_score(results, vina)
        scores['cpd'] = i
        scores['tgt'] = target
        scores.to_csv(scores_path, mode='a', header=None, index=None)

def test():
    import tools
    wt = tools.fasta_to_series('bm3-wt.fasta')[0]
    #prot = protein.protein(pdb_path = '../data/clean/3hf2_clean.pdb',seq = wt)
    df = pd.read_csv('../data/cpds/HRAC_Herbicides.csv')
    vina = Vina(protein='../data/clean/3hf2_clean.pdb')

    vina.read_screen(df)
    vina.dock(df) # should input df here

if __name__ == '__main__':
    test()
