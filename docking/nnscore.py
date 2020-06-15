import oddt
from oddt.scoring.functions import NNScore

# train nn score and pickle

def main():
    path = '../data/clean/1jme_clean.pdb'
    protein = next(oddt.toolkit.readfile('pdb', path))
    protein.addh()

    nn = NNScore.nnscore(protein, n_jobs = 1)



if __name__ == '__main__':
    main()
