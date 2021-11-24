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


pybel.ob.obErrorLog.SetOutputLevel(0)

class pdb_fns:
    def clean_pdb(struc,
            save_path,
            keep = [],
            chain_selection = None):
        if isinstance(keep, str):
            keep = [keep]
        structure = PandasPdb().read_pdb(struc)
        atoms = structure.df['ATOM'].copy()
        hetatms = structure.df['HETATM'].copy()

        uniq_chains = atoms['chain_id'].unique()
        uniq_chains_het = hetatms['chain_id'].unique()
        if chain_selection is None:
            chain_selection = uniq_chains[0]
        if len(uniq_chains_het) == 1:
            chain_selection_het = uniq_chains_het[0]
        else:
            chain_selection_het = chain_selection
        atoms = atoms.loc[atoms['chain_id'] == chain_selection,:]
        hetatms = hetatms.loc[hetatms['chain_id'] == chain_selection_het,:]
        het_garbage = [i for i in hetatms['residue_name'].unique() if i not in keep]
        hetatms = hetatms.loc[hetatms['residue_name'].isin(het_garbage) == False,:]
        structure.df['ATOM'] = atoms
        structure.df['HETATM'] = hetatms
        structure.to_pdb(save_path)
        return save_path

    def get_seq(struc):
        # from vdsl - less error prone
        structure = PandasPdb().read_pdb(struc)
        sequences = structure.amino3to1() # cols = ['chain_id', 'residue_name']
        seqs = [''.join(sequences.loc[sequences['chain_id'] == i,'residue_name'].to_list()) for i in sequences['chain_id'].unique()]
        
        return seqs[0] if len(seqs) == 1 else seqs
        #structure = PandasPdb().read_pdb(struc)
        #sequences = structure.amino3to1() # cols = ['chain_id', 'residue_name']
        #seqs = [''.join(sequences.loc[sequences['chain_id'] == i,'residue_name'].to_list()) for i in sequences['chain_id'].unique()]
        #return seqs[0] if len(seqs) == 1 else seqs

    def draw_box(struc, key_sites):
        receptor = PandasPdb().read_pdb(struc)
        df = receptor.df['ATOM']
        target_site = df.loc[df['residue_number'].isin(key_sites),:]
        coords = target_site.loc[:,['x_coord','y_coord','z_coord']]
        center = coords.mean(axis=0)
        sizes = (coords.max(axis=0) - coords.min(axis=0)) + 5
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
        m.OBMol.StripSalts()
        m.addh()
        m.make3D()
        m.write('pdbqt',save_path, overwrite=True)
        return save_path

    def pdbqt_to_pdb(pdbqt, save_path):
        m = list(pybel.readfile('pdbqt',pdbqt))
        assert len(m) == 1
        m = m[0]
        # already 3d
        m.write('pdb', save_path, overwrite=True)





def aln(s1, s2):
    aln1, aln2 = nw.global_align(s1,s2)
    return aln1, aln2

def diff(s1,s2):
    # s1 - canonical seq
    # s2 - pdb seq
    # offset needed to map aligned positions to pyrosetta positions
    # todo: test where len(s1) < len(s2)
    offset = lambda s, idx : sum([i == '-' for i in s[:idx]])
    return {i - offset(s2,i):{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s2,s1) ,1) if x != y and x != '-' and y != '-'}




