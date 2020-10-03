import os
import pandas as pd
import numpy as np
import multiprocessing
import io
import random
from string import ascii_lowercase

import oddt
from oddt import docking
from oddt import scoring
#from oddt.scoring.functions import NNScore

from biopandas.pdb import PandasPdb

from pyrosetta import pose_from_pdb
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue
import tools
#from enz import NNScore_pdbbind2016.pickle


class protein():
    '''
    object for protein mutation and refolding
    wraps a pyrosetta pose and some cleaning functions

    params: pdb_path, sequence (optional)
    # sequence used for canonical numbering

    example:
    import enz
    seq = <canonical amino acid sequence>
    prot = enz.protein(<pdb path>,seq) # initialise

    prot.mutate(87, 'V') #X87V

    for i in [42, 46, 50, 264, 330]:
        # alanine scan
        prot.mutate(i, 'A')

    prot.refold() # currently: side-chain repacking; planned: loop & flexible region remodelling w/ cyclic coordinate descent

    prot.dump('new.pdb') # save new structure
    '''
    def __init__(self, pdb_path, seq = None):
        self.pdb_path = pdb_path
        self.cache = '__enz-cache__/__protein-cache__' # todo: resolve clashes
        os.makedirs(self.cache, exist_ok=True)
        self._clean_pdb()
        self.pdb_seq = tools.pdb_to_seq(os.path.join(self.cache,'clean.pdb'))
        if seq != None:
            self.seq = seq
        else:
            self.seq = self.pdb_seq
        self.map = tools.map_sequences(self.seq, self.pdb_seq)

    def _clean_pdb(self):
        pdb = tools.clean_pdb(self.pdb_path)
        path_to_clean_pdb = os.path.join(self.cache, 'clean.pdb')
        pdb.to_pdb(path_to_clean_pdb)
        self.path_to_clean_pdb = path_to_clean_pdb # todo: move to tools

    def mutate(self, pos, aa):
        # only changes internal sequence
        # input: self.seq position ; amino acid letter
        seq = list(self.seq) # lists are mutable
        seq[pos] = aa
        seq = ''.join(seq)
        self.seq = seq

    def _mutate_pose(self, pose):
            # mutate at diff between self.seq and self.pdb_seq
            diff = tools.diff(self.pdb_seq,self.seq)
            for i in diff:
                a,b, idx = diff[i]['from'], diff[i]['to'], self.map[i]
                mutate_residue(pose, idx, b, pack_radius = 5.0)

    def refold(self):
        '''
        substitute & repack new sidechains
        '''
        # n decoys?
        # n cpus?
        pyrosetta_init(silent=True)
        if hasattr(self,'path_to_clean_pdb'):
            pose = pose_from_pdb(self.path_to_clean_pdb)
            # mutate
            self._mutate_pose(pose)
            # loop remodel
            self.pose = pose

    def save(self, path):
        if hasattr(self, 'pose'):
            self.pose.dump_pdb(os.path.expanduser(path))
        else:
            self.refold()
            self.pose.dump_pdb(os.path.expanduser(path))
        self._check_dump(os.path.expanduser(path))

    def _check_dump(self,path):
        # check if file ends with END,
        # add END if not
        # otherwise, oddt has trouble reading file
        with open(os.path.expanduser(path),'r') as f:
            file = f.readlines()
        if 'END' not in file[-1]:
            # todo: check other lines for 'END'
            # blank lines not allowed
            file[-1] = file[-1].replace('\n','')
            file.append('END')
            # re-write
            with open(path,'w') as f:
                f.writelines(file) # todo: move to tools

class vina():
    def __init__(self,
     protein = None,
     center = None,
     box_dims = None,
     active_site_aas=None,
     ncpus=None,
     exhaustiveness=8):

        self.cache = '__enz-cache__/__vina-cache__'
        os.makedirs(self.cache, exist_ok=True)
        self.oddt_receptor, self.cache_path = self._process_protein(protein)

        self.box_dims = box_dims
        self.center = center
        self.active_site_aas = active_site_aas
        self.ncpus = ncpus
        self.exhaustiveness = exhaustiveness
        self.vina = self._init_vina()
        self.nn = self.construct_nn()

    def _init_vina(self):
        if self.active_site_aas != None:
            center, box_dims = self._get_centre_pdb()
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

    def _get_centre_pdb(self):
        pdb = PandasPdb().read_pdb(self.cache_path)
        df = pdb.df['ATOM']
        def get_CA(df, aa_num):
            aa = df.loc[df['residue_number']==aa_num,:]
            CA = aa.loc[aa['atom_name'] == 'CA', ['x_coord', 'y_coord', 'z_coord']].values.reshape(-1)
            return CA
        coords = np.array([get_CA(df, i) for i in self.active_site_aas])
        center = coords.mean(axis=0)
        def box_dims(coords, dim):
            return np.abs((coords[:,0].min() - coords[:,0].max()))
        box_dimensions = (box_dims(coords,0), box_dims(coords,1), box_dims(coords,2))
        return center, box_dimensions

    def _process_protein(self,prot):
        # type: pdb path, enz.protein.protein
        # return oddt object
        if isinstance(prot, protein):
            oddt_protein = self._read_protein(prot)
        else:
            # assuming protein is pdb path
            oddt_protein = self._read_pdb(prot)
        path = os.path.join(self.cache, 'receptor.pdb')
        oddt_protein.write('pdb', path, overwrite=True)
        return oddt_protein, path

    def _read_protein(self,protein):
        # if enz.protein: dump in cache
        path = os.path.join(self.cache, 'receptor.pdb')
        protein.save(path) # save as pdb
        oddt_protein = self._read_pdb(path)
        return oddt_protein

    def _read_pdb(self,path):
        oddt_protein = next(oddt.toolkit.readfile('pdb', path))
        oddt_protein.protein = True # register that it's a protein
        oddt_protein.addh()
        return oddt_protein

    def define_active_site(self):
        pass # todo

    def _read_smiles(self,smiles): # todo: move to tools
        mol = oddt.toolkit.readstring('smi',smiles)
        mol.addh()
        mol.make3D()
        return mol

    def _autodock_score(self, results, vina):
        # miving to results obj
        try:
            scores = pd.concat([pd.Series(i.data) for i in vina.score(results)],
            axis=1,
            join='outer').T
            return scores
        except:
            return None

    def construct_nn(self):
        '''
        nn = NNScore.nnscore()
        nn = nn.load(NNScore_pdbbind2016.pickle)
        nn.set_protein(self.oddt_receptor)
        return nn
        '''
        pass

    def nnscore(self, ligands):
        '''
        if type(ligands) == oddt.toolkits.ob.Molecule:
            score  = float(nn.predict_ligand(ligands).data['nnscore'])
        elif type(ligands) == list:
            score = nn.predict(ligands)
        return score
        '''
        pass

    def dock(self, smiles, name=None, ncpu = None, score_fn = None, save = False):
        # set defaults
        if ncpu == None:
            ncpu = multiprocessing.cpu_count() -1
        if score_fn == None:
            score_fn = self._autodock_score
        if name==None:
            name = 'noname'

        # init vina if not done already
        if not hasattr(self,'vina'):
            self.vina = docking.autodock_vina(self.receptor, n_cpu = ncpu) # oddt.docking

        # main bit
        poses = self.vina.dock(self._read_smiles(smiles))
        scores = self._autodock_score(poses, self.vina)
        scores['cpd'] = name
        r = result(poses=poses, receptor=self.oddt_receptor, autodock_score=scores, name=name)

        # todo: get rid of save in vina obj? in favour of save in results obj
        # save
        if save:
            # save scores
            output_dir = os.path.join(self.cache, name)
            os.makedirs(os.path.expanduser(output_dir), exist_ok=True)
            scores.to_csv(os.path.join(output_dir, 'scores.csv'))
            r.save_poses(output_dir)

        return r

class result:
    def __init__(self, poses, receptor, autodock_score, name):
        # todo make coords df
        self.poses = poses
        self.receptor = receptor
        self.autodock_score = autodock_score
        self.name = name # compound
        self.unique_id = 'result-' + ''.join(random.sample(ascii_lowercase,16)) # todo: check unique
        self.cache = os.path.join('__enz-cache__',self.unique_id)
        if not os.path.exists(os.path.join('__enz-cache__',self.unique_id)):
            # why would it exist anyway? get rid?
            os.mkdir(self.cache) # or os.makedirs ??
        self._cache_data(self.cache)
        self.ligand_coords, self.receptor_coords = self._read_cache()

    def _cache_data(self, dirname):
        self.receptor.write('pdb',os.path.join(dirname,'receptor.pdb'))
        for i,mol in enumerate(self.poses):
            filename = f'vina-{self.name}-{i}.pdb' # need to have receptor name too
            path = os.path.join(dirname,filename)
            mol.write('pdb',path,overwrite=True)

    def _read_cache(self):
        # returns atom df
        ligand_coords = {id:PandasPdb().read_pdb(os.path.join(self.cache, i)).df['HETATM'] for id,i in enumerate(os.listdir(self.cache)) if 'vina' in i}
        receptor_obj = PandasPdb().read_pdb(os.path.join(self.cache,'receptor.pdb'))
        receptor_coords = receptor_obj.df['ATOM'].append(receptor_obj.df['HETATM'])

        return ligand_coords, receptor_coords


    def save_autodock_score(self,dirname):
        pass

    def geometric_score(self,prot_atom,lig_atom):
        '''
        EXPERIMENTAL do not keep

        euclidean distance between atoms.
        prot_atom: atom ID in pdb
        lig_atom: atom id in canonical smiles
        '''
        prot_data = PandasPdb().read_pdb(p.path_to_clean_pdb) # move to protein
        # test
        fe  = data.df['HETATM'].loc[data.df['HETATM']['atom_name']=='FE',['x_coord','y_coord','z_coord']].sample() # 2 the same


def test():
    p = protein('1jme.pdb')
    p.mutate(72,'A')
    print('refolding')
    p.refold()
    print('docking')
    v = vina(p)
    r = v.dock('c1cccccc1')
    print(r)

if __name__ == '__main__':
    test()
