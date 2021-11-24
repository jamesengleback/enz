import tempfile
import os
import shutil
from distutils.spawn import find_executable
import subprocess
import re
from itertools import chain

import pandas as pd

import enz
from enz.utils import pdb_fns, obabel_fns


VINA_EXECUTABLE = find_executable('vina')
VINA_SPLIT_EXECUTABLE = find_executable('vina_split')

def dock(receptor_pdb,
         smiles,
         target_sites,
         save_path=None,
         keep=[],
         exhaustiveness=8,
         vina_executable=VINA_EXECUTABLE,
         vina_split_executable=VINA_SPLIT_EXECUTABLE
         ):
    # check there's a box
    #if target_sites == []:
    #    raise Exception('no target_site selected')

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
    args = ['--receptor',
            receptor_pdbqt,
            '--ligand',
            ligand_pdbqt,
            '--out',
            raw_vina_results,
            '--exhaustiveness',
            exhaustiveness]
    
    box = pdb_fns.draw_box(clean_receptor_pdb, target_sites)
    for i,j in zip(box.keys(), box.values()):
        args.append(i)
        args.append(j)
    args = list(map(str, args))
    p1 = subprocess.check_output([vina_executable] + args)
    docking_scores = extract_scores(p1.decode())
    poses = vina_split(raw_vina_results, vina_split_executable)
    results = Results(clean_receptor_pdb, 
                      [os.path.join(poses, i) for i in os.listdir(poses)], 
                      docking_scores)
    
    return results

def vina_split(raw_vina_results, vina_split_executable=VINA_SPLIT_EXECUTABLE):
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

class Pose:
    def __init__(self, 
                 **kwargs
                 ):
        self.__dict__ = {**self.__dict__, **kwargs}
class Results:
    '''
    poses & score df
    '''
    def __init__(self, 
                 receptor, 
                 poses, 
                 vina_scores
                 ):
        self.poses =  {f"mode{get_no(i)}":enz.Mol(i) for i in poses}
        self.receptor = receptor # path to clean pdb
        self.vina_scores = vina_scores.astype(float)
        get_no = lambda path : re.findall('\d+', os.path.basename(path)) # re.findall('\d+',os.path.basename(path))[0]
        self.dictionary = {f"mode{get_no(i)}":{'mol':i, 'affinity':j} \
                for i,j in zip(self.poses.values(), 
                               self.vina_scores['affinity (kcal/mol)'])}
    def save(self, save_path):
        os.makedirs(save_path, exist_ok = True)
        self.vina_scores.to_csv(os.path.join(save_path, 'scores.csv'))
        for i, j in enumerate(self.poses, 1):
            pose_i = self.poses[j]
            pose_i.save(os.path.join(save_path, f'mode{i}.pdb'))
        # saves pdb
        shutil.copyfile(self.receptor, os.path.join(save_path, 'clean_receptor.pdb'))
    def __len__(self):
        return len(self.poses)
    def __getitem__(self, idx):
        return list(self.dictionary.values())[idx] 
    def __repr__(self):
        return f'enz.Results {self.poses}'
