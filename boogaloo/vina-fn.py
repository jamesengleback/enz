import sys
import os
import subprocess
import shutil
import tempfile
from distutils.spawn import find_executable
import re
from openbabel import pybel
import pandas as pd
from biopandas.pdb import PandasPdb

# silence non-critical pybel warnings
pybel.ob.obErrorLog.SetOutputLevel(0)

# executable paths
VINA_EX = '/home/james/src/autodock_vina_1_1_2_linux_x86/bin/vina'
VINA_SPLIT_EX = '/home/james/src/autodock_vina_1_1_2_linux_x86/bin/vina_split' # vina split process the vina output

def dock(receptor, 
        ligand, 
        savepath = None, 
        exhaustiveness=8, 
        cofactors = [], 
        target_residues = [],
        vina_executable = find_executable('vina'),
        vina_split_executable = find_executable('vina_split')):

    # savepath is dir name to dump results in
    CACHE = tempfile.mkdtemp()
    # input prep
    clean_receptor_pdb = clean_pdb(receptor, os.path.join(CACHE, f'{os.path.basename(receptor)}.clean'), cofactors = cofactors)
    receptor_pdbqt = pdb_to_pdbqt(clean_receptor_pdb, os.path.join(CACHE,'receptor.pdbqt'))
    ligand_pdbqt = smiles_to_pdbqt(ligand, os.path.join(CACHE,'ligand.pdbqt'))
    # prep vina arguments 
    outpath = os.path.join(CACHE, 'results.vina')
    args = {'--receptor':receptor_pdbqt,
                '--ligand':ligand_pdbqt,
                '--center_x':10,
                '--center_y':10,
                '--center_z':10,
                '--size_x':10,
                '--size_y':10,
                '--size_z':10,
                '--out':outpath,
                '--exhaustiveness':exhaustiveness}
    # draw vina search box
    if target_residues != []:
        box_dims = get_box_size(clean_receptor_pdb, target_residues)
        for i in box_dims:
            args[i] = box_dims[i]
    # assemble list of args
    args_list_vina = [vina_executable]
    for i in args:
        if args[i] != None:
            args_list_vina.append(i)
            args_list_vina.append(str(args[i]))
    # execute
    p1 = subprocess.check_output(args_list_vina)
    docking_scores = extract_scores(p1.decode())
    # save poses and score csv
    if savepath != None:
        # vina_split
        args_list_vina_split = [vina_split_executable, '--input', outpath]
        p2 = subprocess.Popen(args_list_vina_split, stdout=subprocess.DEVNULL)
        p2.wait()
        # convert pdbqt files 
        docking_poses_pdbqt = [os.path.join(CACHE,i) for i in os.listdir(CACHE) if 'vina_ligand_' in i]
        docking_poses_pdb = [save_pose(i, os.path.join(CACHE, f'pose_{j}.pdb')) for j,i in enumerate(docking_poses_pdbqt)]
        # move to savepath
        os.makedirs(savepath, exist_ok=True)
        # remove existing contents
        for i in os.listdir(savepath):
            os.remove(os.path.join(savepath,i))
        for i in docking_poses_pdb:
            shutil.copyfile(i, os.path.join(savepath, os.path.basename(i))) # copy to inside savepath dir
        shutil.copyfile(clean_receptor_pdb, os.path.join(savepath, os.path.basename(clean_receptor_pdb)))
        docking_scores.to_csv(os.path.join(savepath, 'docking-scores.csv')) # and save scores
    return docking_scores


def extract_scores(text):
    # extract scores from vina output
    text = text.split('\n')
    table_start = ['---+--' in i for i in text].index(True) + 1
    table = []
    for row in text[table_start:]:
        items = row.split()
        is_all_ints = lambda l : sum([re.search('-?\d+', i) is not None for i in l])
        if len(items) == 4 and is_all_ints(items):
            table.append(dict(zip(['mode','affinity (kcal/mol)', 'dist from best mode - rmsd - ub','dist from best mode - lb'], items)))
    return pd.DataFrame(table)

def save_pose(pose,savepath):
    m = list(pybel.readfile('pdbqt',pose))[0]
    m.write('pdb',savepath)
    return savepath

def pdb_to_pdbqt(receptor, savepath):
    m = list(pybel.readfile('pdb',receptor))
    # pybel.readfile returns generator of mols in file
    # hopefully just 1
    assert len(m) == 1
    m = m[0]
    m.addh()
    m.write('pdbqt',savepath, opt={'r':True}, overwrite=True)
    # opt r save pdbqt as rigid, prevents parse errors
    return savepath

def smiles_to_pdbqt(smiles, savepath):
    m = pybel.readstring('smi',smiles)
    m.addh() # redundant - included in make3D()
    m.make3D()
    m.write('pdbqt', savepath)
    return savepath

def clean_pdb(receptor_pdb, savepath, cofactors = [], chain_selection = 'A'):
    # remove: water, duplicate chains
    # assumes that chains are the same, picks one
    receptor = PandasPdb().read_pdb(receptor_pdb)
    receptor.df['ATOM'] = receptor.df['ATOM'].loc[receptor.df['ATOM']['chain_id'] == chain_selection, :]
    throw_away = [i for i in receptor.df['HETATM']['residue_name'].unique() if i not in cofactors]
    receptor.df['HETATM'] = receptor.df['HETATM'].loc[receptor.df['HETATM']['residue_name'].isin(throw_away) == False, :]
    receptor.df['HETATM'] = receptor.df['HETATM'].loc[receptor.df['HETATM']['chain_id'] == chain_selection, :]
    receptor.to_pdb(savepath)
    return savepath


def get_box_size(receptor_pdb, target_residues = []):
    receptor = PandasPdb().read_pdb(receptor_pdb)
    df = receptor.df['ATOM']
    target_site = df.loc[df['residue_number'].isin(target_residues),:]
    coords = target_site.loc[:,['x_coord','y_coord','z_coord']]
    center = coords.mean(axis=0)
    sizes = (coords.max(axis=0) - coords.min(axis=0)) * 1.2
    box_size = {'--center_x':center['x_coord'],
                '--center_y':center['y_coord'],
                '--center_z':center['z_coord'],
                '--size_x':sizes['x_coord'],
                '--size_y':sizes['y_coord'],
                '--size_z':sizes['z_coord']}
    return box_size

def test():
    print(dock('test_data/4key.pdb', 'CCCCCCCCC=O', "results", cofactors='HEM', target_residues=[49,51,82,97,400,300,188,181]))
    # get_box_size('test_data/4key.pdb', [45,87,82,400])
    # clean_pdb('test_data/4key.pdb', CACHE = '.',cofactors=['HEM'])
    # print(get_box_size(clean_pdb('test_data/4key.pdb', cofactors = ['HEM']), target_residues = [49,51,82,97,400,300,188,181]))

if __name__ == '__main__':
    test()
