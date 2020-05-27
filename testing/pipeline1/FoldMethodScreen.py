import pandas as pd
from Bio import pairwise2
import re
from tqdm import tqdm
import os
import random

import pyrosetta
from pyrosetta.rosetta.core.scoring import CA_rmsd, rms_at_all_corresponding_atoms

def renumber(s1,s2, return_both = False):
    alns = pairwise2.align.globalxx(s1,s2,penalize_end_gaps = False,
    one_alignment_only=True)
    s1, s2, _,_,_ = alns[0]
    s1 = re.sub('\w-\w','',s1)
    s2 = re.sub('\w-\w','',s2)
    if return_both:
        return dict(enumerate(s1,0)), dict(enumerate(s2,0))
    else:
        return dict(enumerate(s2,0))

def MapMutations(target,mutant):
    # mutant doesn't have mutations yet
    s1,s2 = target.sequence(), mutant.sequence()
    d1,d2 = renumber(s1,s2, return_both=True)
    mutations = {}
    for i,j in zip(d1,d2):
        if d1[i] != d2[j] and d1[i] != '-' and d2[j] != '-':
            mutations[i] = {'from':d1[i], 'to':d2[j]}

    for i in mutations:
        print('\x1b[31m') #red
        print(mutations[i]['from'], i,mutations[i]['to'])
        pyrosetta.toolbox.mutate_residue(mutant, i,mutations[i]['to'])

def Score(target, mutant):
    # need to also find side chain score
    #atom_rmsd = rms_at_all_corresponding_atoms(target, mutant)
    ca_rmsd = CA_rmsd(target, mutant)
    return ca_rmsd

def PickStructure(path):
    files = [i for i in os.listdir(path) if '.pdb' in i]
    folder = os.path.dirname(path)
    pick = random.sample(files,1)[0]
    return os.path.join(folder, pick)

def PackSideChainsSimple(pose):
    print('\x1b[31m') #red
    print('Packing side chains')
    task_pack = pyrosetta.standard_packer_task(pose)
    task_pack.or_include_current(True)  # considers the original sidechains
    sfxn = pyrosetta.pyrosetta.get_fa_scorefxn()
    pack_mover =  pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sfxn, task_pack)
    pack_mover.apply(pose)

def BackrubMinimise(pose):
    backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
    minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    print('\x1b[31m')
    print('Backrub + minimise')
    for i in tqdm(range(100)):
        backrub.apply(pose)
        minimize.apply(pose)

def main():
    pyrosetta.init()

    # get target and template
    template_file = '../../data/clean/1jme_clean.pdb' # wt???
    target_file = PickStructure('../../data/clean/')
    print('\x1b[31m') #red
    print(f'Target = {os.path.basename(target_file)} \t template = {os.path.basename(template_file)}')
    template = pyrosetta.pose_from_pdb(template_file)
    target = pyrosetta.pose_from_pdb(target_file)

    # make mutant from template
    mutant_pose = pyrosetta.Pose()
    mutant_pose.assign(template)
    MapMutations(target, mutant_pose)

    # score before - rmsd:target
    print('\x1b[31m') #red
    print('score = ',Score(target, mutant_pose))

    # test fold here
    print('\x1b[31m') #red
    BackrubMinimise(mutant_pose)

    # score after
    print('\x1b[31m') #red
    print('score = ',Score(target, mutant_pose))

if __name__ == '__main__':
    main()
