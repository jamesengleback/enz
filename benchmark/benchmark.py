import os 
import re
import random 
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
import enz 
import bm3

def evaluate(template, target, mutations, smiles, ligand_name):
    p_target = enz.protein(target, cofactors=['HEM', ligand_name])
    p = enz.protein(template, cofactors=['HEM'], seq = bm3.BM3_WT)
    if mutations is not np.nan:
        for i in mutations.split(','): # format A82F,etc
            if i != '':
                i = i.strip()
                p.mutate(int(re.findall('([0-9]+)', i)[0]), i[-1])
    else:
        print('mutations', mutations)
    p.refold()
    results = p.dock(smiles, target_residues=bm3.ACTIVE_SITE_AAS)
    # save
    pdb_id = os.path.basename(target).split('.')[0]
    results.save(os.path.join('results', pdb_id))
    p_target.save(os.path.join('results', pdb_id, f'{pdb_id}_.pdb'))


def score(target, results):
    pass

def rmsd(a, b):
    pass

def rotate(coords, x, y, z):
    pass

def align(a, b):
    pass


def test_rotation():
    p = enz.protein('data/1SMJ.pdb')


def main():
    results = os.makedirs('results', exist_ok=True)
    df = pd.read_csv('strucs.tsv', delimiter='\t') 
    combos = random.choices(list(itertools.permutations(df.index, r=2)), k=10)
    for i, j in tqdm(combos):
        template = df.iloc[i,:]
        target = df.iloc[j,:]
        mutations = target['mutations']
        try:
            evaluate(template = template['path'], 
                    target = target['path'], 
                    mutations = mutations, 
                    smiles = target['ligand smiles'],
                    ligand_name = target['ligand id'])
        except Exception as e:
            print(mutations, e)
        


if __name__ == '__main__':
    main()
