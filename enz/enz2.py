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
                chain=None):
        # TODO: select chain
        self.PDB = os.path.abspath(pdb)
        self.ID = uniqueID()
        self.CACHE = tempfile.mkdtemp(prefix='enzp-')
        self.STRUCTURE =  os.path.join(self.CACHE, self.ID + '.pdb') # path
        self._cleanPDB() # save to self.STRUCTURE
        self.seq = self.PDBSEQ if seq == None else seq # self.pdbSeq set in self._cleanPDB
    @property
    def df(self):
        return pd.concat([PandasPdb().read_pdb(self.STRUCTURE).df[i] for i in ['ATOM','HETATM']]).reset_index(drop=True)
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

class vina():
    def __init__(self,
                 pdb,
                 active_site_aas=None,
                 center = None,
                 box_dims = None,
                 ncpus=cpu_count() -1,
                 exhaustiveness=8):
        self.PDB = os.path.abspath(pdb)
        self.ID = uniqueID()
        self.CACHE = tempfile.mkdtemp(prefix='enzv-') #os.path.join('__enz__', self.ID) # tempfile.mkdtemp
        self.STRUCTURE =  os.path.join(self.CACHE, self.ID + '.pdb') # path
        self.active_site_aas = active_site_aas
        self.ncpus = ncpus
        self.exhaustiveness = exhaustiveness

        self._cleanPDB() # problem for heterodimers
        self.vina = self._init_vina()
    @property
    def df(self):
        return pd.concat(
        [PandasPdb().read_pdb(self.j).df[i] for i in ['ATOM','HETATM'] for j in os.listdir(self.CACHE)]).reset_index(drop=True)
    def _cleanPDB(self):
        data = PandasPdb().read_pdb(self.PDB)
        atoms = data.df['ATOM']
        hetatms = data.df['HETATM']
        firstChain = atoms['chain_id'].unique()[0]
        data.df['ATOM'] = atoms.loc[atoms['chain_id'] == firstChain,:]
        data.df['HETATM'] = hetatms.loc[hetatms['residue_name'] != 'HOH',:].loc[hetatms['chain_id'] == firstChain,:]
        data.to_pdb(self.STRUCTURE)
    def _init_vina(self):
        if self.active_site_aas != None:
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
        pass
    def dock(self, smiles, name = None):
        poses = self.vina.dock(self._read_smiles(smiles))
        scores = pd.concat([pd.Series(i.data) for i in self.vina.score(poses,
                        self._read_pdb(self.STRUCTURE))], axis=1, join='outer').T # ValueError: No objects to concatenate
        newIdx = []
        for i,(pose,score) in enumerate(zip(poses, scores['vina_affinity'])):
            savepath = os.path.join(self.CACHE, f'pose-{i}-aff{score}.pdb')
            pose.write('pdb', savepath)
            newIdx.append(savepath)
        scores.index = newIdx
        scores.to_csv(os.path.join(self.CACHE, 'vinaScores.csv'))
        return scores
    def save(self, path):
        shutil.copytree(self.CACHE, path)

def aln(s1, s2):
    aln = global_pairwise_align_protein(skbioProtein(s1),skbioProtein(s2))[0]
    aln_s1 = ''.join([i.decode() for i in aln[0].values])
    aln_s2 = ''.join([i.decode() for i in aln.loc[1].values])
    return aln_s1, aln_s2

def diff(s1,s2):
    return {i:{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s1,s2)) if x != y and x != '-' and y != '-'}
