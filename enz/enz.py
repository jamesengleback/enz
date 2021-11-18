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
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue


PYROSETTA_INIT = False
pybel.ob.obErrorLog.SetOutputLevel(0)

class mol:
    '''
    protein & results poses inherit from this class
    '''
    def __init__(self, struc):
        self.struc = struc
    @property
    def df(self):
        data = PandasPdb().read_pdb(self.struc)
        df = data.df['ATOM'].append(data.df['HETATM']) # all atoms
        return df[['element_symbol', 'residue_number', 'residue_name', 'atom_name', 
            'atom_number', 'x_coord', 'y_coord', 'z_coord']]
    def save(self, save_path):
        shutil.copyfile(self.struc, save_path)


class protein(mol):
    '''
    stores clean protein sequence & structure
    cleaning removes waters and hetatm residues not in keep
    used for mutant structre prediction and small molecule docking
    example:
    >>> import enz
    >>> wt = 'MTIKEM...'
    >>> 
    >>> for i in [75,87,330, 263]:
    >>>    p = enz.protein('enz/data/4key.pdb', seq=wt)
    >>>    p.mutate(i, 'A') # alanine scan
    >>>    r = p.dock('CCCCCCC=O', target_sites=[82,87,330,400,51]) # specify docking site
    >>>    r.save(f'bm3_{p.seq[i]}iA') # e.g. bm3_F87A

    '''
    def __init__(self, struc, seq = None, keep = [], key_sites = [], tmp_suffix=''):
        super().__init__(struc)
        self.key_sites = key_sites
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
        self.struc = folds.fold(self.struc, mutations, pack_radius = 5)
    
    def dock(self, 
            smiles, 
            save_path = None,
            target_sites = None,
            exhaustiveness = 8):
        if target_sites == None:
            target_sites = self.key_sites
        results = vina.dock(self.struc,
                    smiles,
                    save_path = save_path,
                    keep = self.keep,
                    target_sites = target_sites,
                    exhaustiveness = exhaustiveness)
        return results

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


class vina:
    def dock(receptor_pdb,
            smiles,
            save_path = None,
            keep = [],
            target_sites = [],
            exhaustiveness=8,
            vina_executable = find_executable('vina'),
            vina_split_executable = find_executable('vina_split')):
        # check there's a box
        if target_sites == []:
            raise Exception('no target residues selected')

        CACHE = tempfile.mkdtemp(suffix='_enz')
        raw_vina_results = os.path.join(CACHE, 'vina.result')
        # todo : if not clean  - check if dock from protein object
        clean_receptor_pdb = pdb_fns.clean_pdb(receptor_pdb,
                os.path.join(CACHE,
                    f'{os.path.basename(receptor_pdb)}.clean'),
                    keep = keep)
        receptor_pdbqt = obabel_fns.pdb_to_pdbqt(clean_receptor_pdb,
                os.path.join(CACHE,'receptor.pdbqt'))
        ligand_pdbqt = obabel_fns.smiles_to_pdbqt(smiles,
                os.path.join(CACHE,'ligand.pdbqt'))
        args = {'--receptor':receptor_pdbqt,
                    '--ligand':ligand_pdbqt,
                    '--out':raw_vina_results,
                    '--exhaustiveness':exhaustiveness}
        
        args.update(pdb_fns.draw_box(clean_receptor_pdb,
            target_sites)) # add box dims to args

        args_list_vina = [vina_executable] + [str(i) for i in chain.from_iterable(args.items())]

        # execute
        p1 = subprocess.check_output(args_list_vina)
        
        # create results object
        # clean_receptor_pdb to pdb
        docking_scores = vina.extract_scores(p1.decode())
        poses = vina.vina_split(raw_vina_results, vina_split_executable)
        results = vina.results(clean_receptor_pdb, [os.path.join(poses, i) for i in os.listdir(poses)], docking_scores)
        
        return results

    def vina_split(raw_vina_results, vina_split_executable):
        # vina_split
        args_list_vina_split = [vina_split_executable, '--input', raw_vina_results]
        p = subprocess.Popen(args_list_vina_split, stdout=subprocess.DEVNULL)
        p.wait() # outputs odbqt files
        results_dir = os.path.dirname(raw_vina_results)
        poses = [os.path.join(results_dir, i) for i in os.listdir(results_dir) if 'vina.result_ligand' in i]
        clean_results = os.path.join(results_dir, 'pose_pdbs')
        os.makedirs(clean_results, exist_ok=True)
        for i in poses:
            save_path = os.path.basename(i).replace('pdbqt','pdb')
            obabel_fns.pdbqt_to_pdb(i, os.path.join(clean_results, save_path))
        return clean_results # path to foder containing pdb of poses

    def extract_scores(text):
        # extract scores from vina output
        text = text.split('\n')
        table_start = ['---+--' in i for i in text].index(True) + 1
        is_all_ints = lambda l : sum([re.search('-?\d+', i) is not None for i in l])
        table = []
        for row in text[table_start:]:
            items = row.split()
            if len(items) == 4 and is_all_ints(items):
                table.append(dict(zip(['mode','affinity (kcal/mol)', 'dist from best mode - rmsd - ub','dist from best mode - lb'], items)))
        return pd.DataFrame(table)

    class results:
        '''
        poses & score df
        '''
        def __init__(self, receptor, poses, scores):
            self.poses = {int(re.findall('\d+',os.path.basename(i))[0])\
                    :mol(i) for i in poses}
            self.receptor = receptor # path to clean pdb
            self.scores = scores.astype(float)
            self.dictionary = {os.path.basename(i.struc):{'mol':i, 'affinity':j} for i,j in zip(self.poses.values(), self.scores['affinity (kcal/mol)'])}
        def save(self, save_path):
            os.makedirs(save_path, exist_ok = True)
            self.scores.to_csv(os.path.join(save_path, 'scores.csv'))
            for i, j in enumerate(self.poses, 1):
                pose_i = self.poses[j]
                pose_i.save(os.path.join(save_path, f'mode{i}.pdb'))
            # saves pdb
            shutil.copyfile(self.receptor, os.path.join(save_path, 'clean_receptor.pdb'))



class utils:
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



PYROSETTA_INIT = False

class folds:
    def init():
        global PYROSETTA_INIT
        if  not PYROSETTA_INIT:
            pyrosetta_init(silent=True)
            PYROSETTA_INIT = True

    def getpose(pdb):
        return pose_from_pdb(pdb)

    def savepose(pose):
        tmp = tempfile.mktemp('_enz')
        pose.dump_file(tmp)
        return tmp

    def fold(pdb, mutation_dict, pack_radius = 5):
        folds.init() 
        pose = folds.getpose(pdb)
        pose = folds.fold_repack_mutate(pose, mutation_dict, pack_radius)
        return folds.savepose(pose) # returns tempfile path

    def fold_repack_mutate(pose, mutation_dict, pack_radius = 5):

        for i in mutation_dict:
            print(i, mutation_dict[i])
            mutate_residue(pose, i,
                            mutation_dict[i]['to'].upper(),
                            pack_radius = float(pack_radius))
        return pose

    def fold_ccd(pose, mutation_dict):
        folds.init()
        def detect_loops(pose):
            # phi psi
            phipsi = [[pose.phi(i), pose.psi(i)] for i, _ in enumerate(pose.sequence(), 1)]
            print(phipsi)
        def ccd_loop(pos):
            pass
        loops = detect_loops(pose)
        loops = [] # delet
        for i in mutation_dict:
            if i in loops:
                ccd_loop(i)
        return pose


