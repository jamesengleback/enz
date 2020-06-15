import oddt
from oddt import docking
from oddt.scoring.functions import NNScore
import numpy as np
import pandas as pd
import os
from tqdm import  tqdm

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

def dock_multi(protein, smiles, names):
    # protein must be oddt type
    
    outdir = 'DockingResults'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    scores_path = os.path.join(outdir,'autodock_scores.csv')

    df = pd.DataFrame([],columns=['vina_affinity', 'vina_rmsd_lb', 'vina_rmsd_ub', 'vina_rmsd_input',
       'vina_rmsd_input_min', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion',
       'vina_hydrophobic', 'vina_hydrogen', 'SMILES Atom Order', 'name'])

    df.to_csv(scores_path)

    vina = docking.autodock_vina(protein)

    for i,smi in tqdm(zip(names, smiles), total=len(smiles)):
        mol = read_smiles(smi)
        results = vina.dock(mol)
        save_results(results, i, parent_dir=outdir)
        scores = autodock_score(results, vina)
        scores['name'] = i
        scores.to_csv(scores_path, mode='a')

def test():
    protein = make_protein('../data/clean/1jme_clean.pdb')
    #mol = read_smiles('c1ccccc1')
    #s = dock_single(protein,mol)
    #print(s)
    df = pd.read_csv('../data/cpds/HRAC_Herbicides.csv')
    names = df['Name']
    smiles = df['SMILES']
    dock_multi(protein, smiles,names)


if __name__ == '__main__':
    test()
