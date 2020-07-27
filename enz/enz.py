import os
import pandas as pd
import numpy as np
from tqdm import  tqdm
import multiprocessing
import io

import oddt
from oddt import docking
from oddt import scoring
from oddt.scoring.functions import NNScore

from biopandas.pdb import PandasPdb

import pyrosetta

from enz import tools # doesnt work in testiing
from enz import NNScore_pdbbind2016.pickle


class Protein():
    def __init__(self, pdb_path, seq = None):
        self.pdb_path = pdb_path
        self.cache = '__protein-cache__' # todo: resolve clashes
        os.makedirs(self.cache, exist_ok=True)
        self.clean_pdb()
        self.pdb_seq = tools.pdb_to_seq(os.path.join(self.cache,'clean.pdb'))
        self.seq = seq if seq != None else tools.pdb_to_seq(self.pdb_path)
        self.map = tools.map_sequences(self.seq, self.pdb_seq)

    def clean_pdb(self):
        pdb = tools.clean_pdb(self.pdb_path)
        path_to_clean_pdb = os.path.join(self.cache, 'clean.pdb')
        pdb.to_pdb(path_to_clean_pdb)
        self.path_to_clean_pdb = path_to_clean_pdb # todo: move to tools

    def mutate_seq(self, pos, aa):
        # only changes internal sequence
        # input: self.seq position ; amino acid letter
        seq = list(self.seq) # lists are mutable
        seq[pos] = aa
        seq = ''.join(seq)
        self.seq = seq

    def mutate_pose(self, pose):
            # mutate at diff between self.seq and self.pdb_seq
            diff = tools.diff(self.pdb_seq,self.seq)
            for i in diff:
                a,b, idx = diff[i]['from'], diff[i]['to'], self.map[i]
                pyrosetta.toolbox.mutate_residue(pose, idx, b, pack_radius = 5.0)

    def refold(self):
        '''
        substitute & repack new sidechains
        '''
        # n decoys?
        # n cpus?
        pyrosetta.init(silent=True)
        if hasattr(self,'path_to_clean_pdb'):
            pose = pyrosetta.pose_from_pdb(self.path_to_clean_pdb)
            # mutate
            self.mutate_pose(pose)
            # loop remodel
            self.pose = pose

    def dump(self, path):
        if hasattr(self, 'pose'):
            self.pose.dump_pdb(path)
        else:
            self.refold()
            self.pose.dump_pdb(path)
        self.check_dump(path)

    def check_dump(self,path):
        # check if file ends with END,
        # add END if not
        # otherwise, oddt has trouble reading file
        with open(path,'r') as f:
            file = f.readlines()
        if 'END' not in file[-1]:
            # todo: check other lines for 'END'
            # blank lines not allowed
            file[-1] = file[-1].replace('\n','')
            file.append('END')
            # re-write
            with open(path,'w') as f:
                f.writelines(file) # todo: move to tools

class Vina():
    def __init__(self,
     protein = None,
     center = None,
     box_dims = None,
     acitve_site_aas=None,
     ncpus=None,
     exhaustiveness=8):

        self.cache = '__vina-cache__'
        os.makedirs(self.cache, exist_ok=True)
        self.oddt_receptor, self.cache_path = self.process_protein(protein)

        self.box_dims = box_dims
        self.center = center
        self.acitve_site_aas = acitve_site_aas
        self.ncpus = ncpus
        self.exhaustiveness = exhaustiveness
        self.vina = self.init_vina()
        self.nn = self.construct_nn()

    def init_vina(self):
        if self.acitve_site_aas != None:
            center, box_dims = self.get_centre_pdb()
        else:
            center, box_dims = (0,0,0), (20,20,20)
        if self.ncpus == None:
            self.ncpus = multiprocessing.cpu_count() -1
        if self.exhaustiveness == None:
            self.exhaustiveness = 8

        vina = docking.autodock_vina(self.oddt_receptor,
        n_cpu = self.ncpus,
        exhaustiveness=self.exhaustiveness,
        size=box_dims,
        center=center)
        return vina

    def get_centre_pdb(self):
        pdb = PandasPdb().read_pdb(self.cache_path)
        df = pdb.df['ATOM']
        def get_CA(df, aa_num):
            aa = df.loc[df['residue_number']==aa_num,:]
            CA = aa.loc[aa['atom_name'] == 'CA', ['x_coord', 'y_coord', 'z_coord']].values.reshape(-1)
            return CA
        coords = np.array([get_CA(df, i) for i in self.acitve_site_aas])
        center = coords.mean(axis=0)
        def box_dims(coords, dim):
            return np.abs((coords[:,0].min() - coords[:,0].max()))
        box_dimensions = (box_dims(coords,0), box_dims(coords,1), box_dims(coords,2))
        return center, box_dimensions

    def process_protein(self,prot):
        # type: pdb path, enz.protein.protein
        # return oddt object
        if isinstance(prot, Protein):
            oddt_protein = self.read_protein(prot)
        else:
            # assuming protein is pdb path
            oddt_protein = self.read_pdb(prot)
        path = os.path.join(self.cache, 'receptor.pdb')
        oddt_protein.write('pdb', path, overwrite=True)
        return oddt_protein, path

    def read_protein(self,protein):
        # if enz.protein: dump in cache
        path = os.path.join(self.cache, 'receptor.pdb')
        protein.dump(path) # save as pdb
        oddt_protein = self.read_pdb(path)
        return oddt_protein

    def read_pdb(self,path):
        oddt_protein = next(oddt.toolkit.readfile('pdb', path))
        oddt_protein.protein = True # register that it's a protein
        oddt_protein.addh()
        return oddt_protein

    def define_active_site(self):
        pass # todo

    def read_smiles(self,smiles): # todo: move to tools
        mol = oddt.toolkit.readstring('smi',smiles)
        mol.addh()
        mol.make3D()
        return mol

    def autodock_score(self, results, vina):
        try:
            scores = pd.concat([pd.Series(i.data) for i in vina.score(results)],
            axis=1,
            join='outer').T
            return scores
        except:
            raise EnzError('Auto dock score failed')
            return None

    def construct_nn(self):
        nn = NNScore.nnscore()
        nn = nn.load(NNScore_pdbbind2016.pickle)
        nn.set_protein(self.oddt_receptor)
        return nn

    def nnscore(self, ligands):
        if type(ligands) == oddt.toolkits.ob.Molecule:
            score  = float(nn.predict_ligand(ligands).data['nnscore'])
        elif type(ligands) == list:
            score = nn.predict(ligands)
        return score

    def dock(self, smiles, name=None, ncpu = None, score_fn = None, save = True):
        # set defaults
        if ncpu == None:
            ncpu = multiprocessing.cpu_count() -1
        if score_fn == None:
            score_fn = self.autodock_score
        if name==None:
            name = 'noname'

        # init vina if not done already
        if not hasattr(self,'vina'):
            self.vina = docking.autodock_vina(self.receptor, n_cpu = ncpu) # oddt.docking

        # main bit
        try:
            results = self.vina.dock(self.read_smiles(smiles))
            scores = self.autodock_score(results, self.vina)
            scores['cpd'] = name

            # save
            if save:
                # save scores
                output_dir = os.path.join(self.cache, name)
                os.makedirs(output_dir, exist_ok=True)
                scores.to_csv(os.path.join(output_dir, 'scores.csv'))

                # save docking poses
                for i,mol in enumerate(results):
                    filename = f'vina-{name}-{i}.pdb' # need to have receptor name too
                    path = os.path.join(output_dir,filename)
                    mol.write('pdb',path,overwrite=True)

            return scores, results
        except:
            raise EnzError('Uh oh - error in docking - could be heavy atoms?')
            return None, None


class EnzError(Exception):
    pass


def test():
    bm3_wt = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRF\
    DKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDM\
    TRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADR\
    KASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDP\
    VPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERF\
    ENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGG\
    IPSPSTEQSAKKVRKKGC*'.replace(' ','')


    vina = Vina('../test/data/clean/3ben_clean.pdb',
    acitve_site_aas = [47, 50, 51, 72, 75, 78, 82, 87, 181, 188, 263, 268, 328, 330])


    scores, results  = vina.dock('cccccccc')
    print(scores)
    for i,  mol in enumerate(results):
        mol.write('pdb', f'__{i}.pdb', overwrite=True)

if __name__ == '__main__':
    test()
