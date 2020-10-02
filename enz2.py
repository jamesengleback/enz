import sys
import os
from string import ascii_lowercase, ascii_uppercase
import random
import pandas as pd
from biopandas.pdb import PandasPdb

from pyrosetta import pose_from_pdb
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue

from skbio.alignment import global_pairwise_align_protein

CACHE = '__enz-cache__'
PYROSETTA_INIT = False
os.makedirs(CACHE, exist_ok=True)


class protein:
    '''
    '''
    def __init__(self, pdb = None, seq = None):
        self.PDB = os.path.abspath(pdb)
        self.ID = ''.join(random.sample(ascii_lowercase + ascii_uppercase, 20))
        self.STRUCTURE =  os.path.join(CACHE, self.ID + '.pdb') # path
        self._cleanPDB() # save to self.STRUCTURE
        self.seq = self.PDBSEQ if seq == None else seq # self.pdbSeq set in self._cleanPDB
    @property
    def df(self):
        return pd.concat([PandasPdb().read_pdb(self.STRUCTURE).df[i] for i in ['ATOM','HETATM']]).reset_index(drop=True)
    def _cleanPDB(self):
        data = PandasPdb().read_pdb(self.PDB)
        atoms = data.df['ATOM']
        hetatms = data.df['HETATM']
        data.df['ATOM'] = atoms.loc[atoms['chain_id'] == atoms['chain_id'].unique()[0],:]
        data.df['HETATM'] = hetatms.loc[hetatms['residue_name'] != 'HOH',:]
        data.to_pdb(self.STRUCTURE)
        self.PDBSEQ = ''.join([i for i in data.amino3to1()['residue_name']]) # PDBSEQ
    def _align(self):
        pass
    def mutate(self, idx, aa):
        seq = list(self.seq)
        seq[idx] = aa
        self.seq = ''.join(seq)
    def refold(self):
        pass
    def save(self, path):
        pass


#        if not PYROSETTA_INIT:
#            pyrosetta_init(silent=True)
