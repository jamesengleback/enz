import oddt
from oddt import docking
import numpy as np

def read_smiles(smiles):
    mol = oddt.toolkit.readstring('smi',smiles)
    mol.addh()
    mol.make3D()
    return mol

def find_heme(protein):
    atoms = [i.atomicnum for i in list(protein.atoms)]
    coords = protein.coords
    d = {idx:{'atom':a, 'coords':c} for idx, a, c in zip(range(len(atoms)), atoms, coords)}
    heme_fe = [d[i] for i in d if d[i]['atom'] == 26]
    return next(iter(heme_fe)) # gets it out of a list

def euclidian_distance(a,b):
    return np.linalg.norm(a-b)

def distance_to(mol,target):
    d = {}
    for atom in mol:
        a = np.array(atom.coords)
        d[atom.idx0] = euclidian_distance(a, target['coords'])
        # extract more info
    print(d)


def main():
    path = '../data/clean/1jme_clean.pdb'
    protein = next(oddt.toolkit.readfile('pdb', path))
    protein.addh()

    mol = read_smiles('c1ccccc1')

    vina = docking.autodock_vina(protein)
    results = vina.dock(mol)
    heme_fe = find_heme(protein)
    print(distance_to(results[0], heme_fe))


if __name__ == '__main__':
    main()
