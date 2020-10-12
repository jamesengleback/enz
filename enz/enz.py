import os
import shutil
import tempfile
from string import ascii_lowercase, ascii_uppercase
import random
from multiprocessing import cpu_count

import pandas as pd
from biopandas.pdb import PandasPdb

from pyrosetta import pose_from_pdb
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue

from skbio.alignment import global_pairwise_align_protein
from  skbio import Protein as skbioProtein

import oddt
from oddt import docking

PYROSETTA_INIT = False

def uniqueID():
    return ''.join(random.sample(ascii_lowercase + ascii_uppercase, 20))
class protein:
    '''
    '''
    def __init__(self,
                pdb = None,
                seq = None,
                key_sites = {},
                chain=None):
        # TODO: select chain
        self.PDB = os.path.abspath(pdb)
        self.ID = uniqueID()
        self.CACHE = tempfile.mkdtemp(prefix='enzp-')
        self.STRUCTURE =  os.path.join(self.CACHE, self.ID + '.pdb') # path
        self._cleanPDB() # save to self.STRUCTURE
        self.seq = self.PDBSEQ if seq == None else seq # self.pdbSeq set in self._cleanPDB
        self.KEY_SITES = {i:self.seq[i] for i in key_sites}
    @property
    def df(self):
        return pd.concat([PandasPdb().read_pdb(self.STRUCTURE).df[i] for i in ['ATOM','HETATM']]).reset_index(drop=True)
    @property
    def KEY_SITES_DICT(self):
        if self.KEY_SITES != None:
            return {i:self.seq[i] for i in self.KEY_SITES}
    @property
    def docking_results(self):
        if hasattr(self, 'vina'):
            return self.vina.CACHE
    def _cleanPDB(self):
        data = PandasPdb().read_pdb(self.PDB)
        atoms = data.df['ATOM']
        hetatms = data.df['HETATM']
        firstChain = atoms['chain_id'].unique()[0]
        data.df['ATOM'] = atoms.loc[atoms['chain_id'] == firstChain,:]
        data.df['HETATM'] = hetatms.loc[hetatms['residue_name'] != 'HOH',:].loc[hetatms['chain_id'] == firstChain,:]
        data.to_pdb(self.STRUCTURE)
        self.PDBSEQ = ''.join([i for i in data.amino3to1()['residue_name']]) # PDBSEQ
    def mutate(self, idx, aa):
        seq = list(self.seq)
        seq[idx] = aa
        self.seq = ''.join(seq)
    def refold(self):
        aln1, aln2 = aln(self.seq, self.PDBSEQ)
        mutations = diff(aln1, aln2)
        if  not PYROSETTA_INIT:
            pyrosetta_init(silent=True)
            PYROSETTA_INIT = True
        self.pose = pose_from_pdb(self.STRUCTURE)
        for i in mutations:
            mutate_residue(self.pose, i, mutations[i]['to'], pack_radius = 5.0)
        self.pose.dump_pdb(self.STRUCTURE)
    def save(self, path):
        shutil.copy(self.STRUCTURE, path)
    def dock(self,smiles, name=None):
        if not hasattr(self, 'vina'):
            self.vina = vina(pdb = self.STRUCTURE,
                            key_sites = self.KEY_SITES,)
        scores = self.vina.dock(smiles,name)
        return scores
    def save_docking_results(self,path):
        if hasattr(self,'vina'):
            self.vina.save(path)



class vina():
    def __init__(self,
                 pdb,
                 key_sites=None,
                 center = None,
                 box_dims = None,
                 ncpus=cpu_count() -1,
                 exhaustiveness=8):
        self.PDB = os.path.abspath(pdb)
        self.ID = uniqueID()
        self.CACHE = tempfile.mkdtemp(prefix='enzv-') #os.path.join('__enz__', self.ID) # tempfile.mkdtemp
        self.STRUCTURE =  os.path.join(self.CACHE, self.ID + '.pdb') # path
        self.KEY_SITES = key_sites
        self.ncpus = ncpus
        self.exhaustiveness = exhaustiveness
        self._cleanPDB() # problem for heterodimers
        self.vina = self._init_vina()
    @property
    def df(self):
        cached_files = [os.path.join(parent,file) for parent, _, files in os.walk(self.CACHE) for file in files]
        return pd.concat([PandasPdb().read_pdb(j).df[i] for i in ['ATOM','HETATM'] for j in cached_files])
    def _cleanPDB(self): # slow?
        data = PandasPdb().read_pdb(self.PDB)
        atoms = data.df['ATOM']
        hetatms = data.df['HETATM']
        firstChain = atoms['chain_id'].unique()[0]
        data.df['ATOM'] = atoms.loc[atoms['chain_id'] == firstChain,:]
        data.df['HETATM'] = hetatms.loc[hetatms['residue_name'] != 'HOH',:].loc[hetatms['chain_id'] == firstChain,:]
        data.to_pdb(self.STRUCTURE)
    def _init_vina(self):
        if self.KEY_SITES != None:
            center, box_dims = self._boxDims()
        else:
            center, box_dims = (0,0,0), (20,20,20)
        vina = docking.autodock_vina(self._read_pdb(self.STRUCTURE),
                                                n_cpu = self.ncpus,
                                                exhaustiveness=self.exhaustiveness,
                                                size=box_dims,
                                                center=center)
        return vina
    def _read_pdb(self,path):
        oddt_protein = next(oddt.toolkit.readfile('pdb', path))
        oddt_protein.protein = True # register that it's a protein
        oddt_protein.addh()
        return oddt_protein
    def _read_smiles(self,smiles): # todo: move to tools
        mol = oddt.toolkit.readstring('smi',smiles)
        mol.addh()
        mol.make3D()
        return mol
    def _boxDims(self):
        df = self.df
        key_sites_coords = df.loc[df['residue_number'].isin(self.KEY_SITES), ['x_coord','y_coord','z_coord']]
        center = key_sites_coords.mean(axis=0)
        x,y,z = [(key_sites_coords[i].max() - key_sites_coords[i].min()) * 1.2 for i in key_sites_coords]
        return center, (x,y,z)

    def dock(self, smiles, name = None):
        poses = self.vina.dock(self._read_smiles(smiles))
        scores = pd.concat([pd.Series(i.data) for i in self.vina.score(poses,
                        self._read_pdb(self.STRUCTURE))], axis=1, join='outer').T # ValueError: No objects to concatenate
        newIdx = []
        for i,(pose,score) in enumerate(zip(poses, scores['vina_affinity'])):
            savepath = os.path.join(self.CACHE, f'pose-{uniqueID()}:aff{score}.pdb')
            pose.write('pdb', savepath)
            newIdx.append(savepath)
        scores['path'] = newIdx
        scores.to_csv(os.path.join(self.CACHE, 'vinaScores.csv'))
        return scores
    def save(self, path):
        shutil.copytree(self.CACHE, path, dirs_exist_ok=True)

def aln(s1, s2):
    aln = global_pairwise_align_protein(skbioProtein(s1),skbioProtein(s2))[0]
    aln_s1 = ''.join([i.decode() for i in aln[0].values])
    aln_s2 = ''.join([i.decode() for i in aln.loc[1].values])
    return aln_s1, aln_s2

def diff(s1,s2):
    return {i:{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s1,s2)) if x != y and x != '-' and y != '-'}
