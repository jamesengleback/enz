import sys
import os
import tempfile

from pyrosetta import pose_from_pdb
from pyrosetta import logger as pyrosetta_logger
from pyrosetta import logging as pyrosetta_logging
from pyrosetta import init as pyrosetta_init
from pyrosetta.toolbox import mutate_residue

import logging
#logging.getLogger("pyrosetta.rosetta").setLevel(0)

PYROSETTA_INIT = False

def init():
    global PYROSETTA_INIT
    if not PYROSETTA_INIT:
        pyrosetta_init(silent=True) # loads of text output
        # v causes threading issues????
        #pyrosetta_init(silent=True, set_logging_handler='logging') # actually silent
        PYROSETTA_INIT = True

def fold(pdb, 
         mutation_dict, 
         pack_radius = 5):
    '''
    main method, combines functions
    '''
    sys.stderr = open(os.devnull, 'w')
    init() 
    pose = getpose(pdb)
    pose = fold_repack_mutate(pose, mutation_dict, pack_radius)
    sys.stdout = sys.__stdout__
    return savepose(pose) # returns tempfile path

def getpose(pdb):
    return pose_from_pdb(pdb)

def savepose(pose):
    tmp = tempfile.mktemp('_enz')
    pose.dump_file(tmp)
    return tmp

def fold_repack_mutate(pose, mutation_dict, pack_radius = 5):

    for i in mutation_dict:
        mutate_residue(pose, i,
                        mutation_dict[i]['to'].upper(),
                        pack_radius = float(pack_radius))
    return pose

def fold_ccd(pose, mutation_dict):
    #init()
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
