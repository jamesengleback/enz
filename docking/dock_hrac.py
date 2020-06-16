from docking_tools import dock_multi, make_protein
import pandas as pd
import os
from tqdm import  tqdm

def main():
    paths = ['../data/clean/3ben_clean.pdb',
    '../data/clean/1smi_clean.pdb',
    '../data/clean/1p0x_clean.pdb',
    '../data/clean/1jme_clean.pdb']

    df = pd.read_csv('../data/cpds/HRAC_Herbicides.csv').sample(1)
    names = df['Name']
    smiles = df['SMILES']

    for i in paths:
        protein = make_protein(i)
        target = os.path.basename(i).split('_')[0]
        print(target)
        dock_multi(protein, smiles,names, target)



if __name__ == '__main__':
    main()
