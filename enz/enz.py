import tempfile
import os
import shutil
from distutils.spawn import find_executable
import subprocess
import re
from itertools import chain

import pandas as pd

from biopandas.pdb import PandasPdb
import nwalign3 as nw
from openbabel import pybel

from pyrosetta import pose_from_pdb
from pyrosetta import logger as pyrosetta_logger
from pyrosetta import logging as pyrosetta_logging
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue

import enz.utils as utils
from enz.utils import pdb_fns, obabel_fns, aln, diff
import enz.vina as vina
import enz.fold as fold

PYROSETTA_INIT = False
pyrosetta_logger.setLevel(0)
pyrosetta_logging.disable()
pybel.ob.obErrorLog.SetOutputLevel(0)

class Mol:
    '''
    protein & results poses inherit from this class
    '''
    def __init__(self, 
                 struc):
        self.struc = struc # pdb path
    @property
    def df(self):
        data = PandasPdb().read_pdb(self.struc)
        df = data.df['ATOM'].append(data.df['HETATM']) # all atoms
        return df[['element_symbol', 'residue_number', 'residue_name', 'atom_name', 
            'atom_number', 'x_coord', 'y_coord', 'z_coord']]
    def save(self, save_path):
        shutil.copyfile(self.struc, save_path)
    def __repr__(self):
        return f'enz.Mol {self.struc}'


class Protein(Mol):
    '''
    stores clean protein sequence & structure
    cleaning removes waters and hetatm residues not in keep
    used for mutant structre prediction and small molecule docking
    example:
    >>> import enz
    >>> wt = 'MTIKEM...'
    >>> 
    >>> for i in [75,87,330, 263]:
    >>>    p = enz.Protein('enz/data/4key.pdb', seq=wt)
    >>>    p.mutate(i, 'A') # alanine scan
    >>>    p.refold()
    >>>    r = p.dock('CCCCCCC=O', target_sites=[82,87,330,400,51]) # specify docking site
    >>>    r.save(f'bm3_{p.seq[i]}iA') # e.g. bm3_F87A

    '''
    def __init__(self, 
                 struc, 
                 seq = None, 
                 keep = [], 
                 tmp_suffix=''):
        super().__init__(struc)
        self.keep = keep
        self.CACHE = tempfile.mkdtemp(suffix=f'{tmp_suffix}_enz')
        self.struc = pdb_fns.clean_pdb(struc = self.struc,
                                      save_path = os.path.join(self.CACHE, 'clean.pdb'),
                                      keep = self.keep)
        self.pdb_seq = pdb_fns.get_seq(self.struc)
        self.seq = self.pdb_seq if seq == None else seq

    def mutate(self, position, aa):
        seq = list(self.seq)
        seq[position] = aa
        self.seq = ''.join(seq)
    
    def refold(self, pack_radius = 5.0):
        aln1, aln2 = utils.aln(self.seq, self.pdb_seq)
        mutations = utils.diff(aln1, aln2)
        self.struc = fold.fold(self.struc, mutations, pack_radius = 5)
    
    def dock(self, 
            smiles, 
            save_path = None,
            target_sites = None,
            exhaustiveness = 8):
        return vina.dock(self.struc,
                         smiles,
                         save_path=save_path,
                         keep=self.keep,
                         target_sites=target_sites,
                         exhaustiveness=exhaustiveness)


