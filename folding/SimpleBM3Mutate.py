import pyrosetta
from pyrosetta import rosetta
from pyrosetta.rosetta.core.scoring import CA_rmsd

from tqdm import tqdm
import argparse
import sys

colors = {'Black': '\x1b[30m',
 'Red': '\x1b[31m',
 'Green': '\x1b[32m',
 'Yellow': '\x1b[33m',
 'Blue': '\x1b[34m',
 'Magenta': '\x1b[35m',
 'Cyan': '\x1b[36m',
 'White': '\x1b[37m',
 'Reset': '\x1b[0m',
 'BrightBlack': '\x1b[30;1m',
 'BrightRed': '\x1b[31;1m',
 'BrightGreen': '\x1b[32;1m',
 'BrightYellow': '\x1b[33;1m',
 'BrightBlue': '\x1b[34;1m',
 'BrightMagenta': '\x1b[35;1m',
 'BrightCyan': '\x1b[36;1m',
 'Bold': '\x1b[1m',
 'Underline': '\x1b[4m',
 'Reversed': '\x1b[7m'}

def ParseMutation(m):
    number = int(''.join([i for i in m if i.isdigit()]))
    letter = [i for i in m if i.isalpha()][-1] # in case they give 2 letters
    # If it's 1BU7, then numbering = normal numbering -2
    return number, letter.upper()

def Classic_relax(pose):
    # classic relax protocol
    #P. Bradley, K. M. S. Misura & D. Baker, “Toward high-resolution de novo structure
    #prediction for small proteins,” Science 309, 1868-1871 (2005)
    # classic relax alternates an energy minimizer and a side chain packer
    ######### Doesn't work very well!!
    sfxn = pyrosetta.ScoreFunction() # default score function
    relax = pyrosetta.rosetta.protocols.relax.ClassicRelax()
    relax.set_scorefxn(sfxn)
    relax.apply(pose)
    return pose

def BackrubMinimise(pose):
    # works ok
    # how many iterations?
    # what else?
    backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
    minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    print(colors['BrightBlue'])
    print('Backrub + minimise')
    for i in tqdm(range(100)):
        backrub.apply(pose)
        minimize.apply(pose)
    return pose

def PackSideChainsSimple(pose):
    print(colors['BrightGreen'])
    print('Packing side chains')
    task_pack = pyrosetta.standard_packer_task(pose)
    task_pack.or_include_current(True)  # considers the original sidechains
    sfxn = pyrosetta.pyrosetta.get_fa_scorefxn()
    pack_mover =  pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sfxn, task_pack)
    pack_mover.apply(pose)
    return pose

def PackSideChains_Complicated(pose, kT = 0.8, cycles = 10):
    starting_pose = Pose()
    starting_pose.assign(pose)

    print(colors['BrightGreen'])
    print('Packing side chains')

    task_pack = pyrosetta.standard_packer_task(pose)
    task_pack.or_include_current(True)  # considers the original sidechains
    sfxn = pyrosetta.pyrosetta.get_fa_scorefxn()
    pack_mover =  pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sfxn, task_pack)
    mc = MonteCarlo(pose, sfxn, kT)
    trial = TrialMover(pack_mover, mc)
    job = RepeatMover(trial, cycles)

    jd = PyJobDistributor(job_output, jobs, sfxn)
    jd.native_pose = starting_pose
    while not jd.job_complete:
        pose.assign(starting_pose)
        mc.reset(pose)
        job.apply(pose)
        mc.recover_low(pose)
        jd.output_decoy(pose)
    return pose

def main(args):

    if args.template != None:
        pyrosetta.init()
        template = pyrosetta.pose_from_pdb(args.template)
    else:
        print(colors['BrightRed'])
        print('No input template specified')
        print('Exiting')
        sys.exit()

    print(colors['BrightYellow'])
    print('Template:\tOk')
    print(colors['Reset'])


    if args.mutation != None:
        mutation_position, mutate_to = ParseMutation(args.mutation)
        print(colors['BrightYellow'])
        print(f'Making Mutant {template.sequence()[mutation_position-1]}{mutation_position}{mutate_to}')
        print(colors['Reset'])
    else:
        print(colors['BrightRed'])
        print('No mutation!')
        print('Exiting')
        sys.exit()


    # mutation stuff
    # copy template
    mutant_pose = pyrosetta.Pose()
    mutant_pose.assign(template)
    mutant_pose.pdb_info().name(f'{template.sequence()[mutation_position-1]}\
    {mutation_position}{mutate_to}')
    #make mutation
    pyrosetta.toolbox.mutate_residue(mutant_pose, mutation_position-1,mutate_to)

    # minimise energy
    mutant_pose = PackSideChains(mutant_pose)
    #dump output
    mutant_pose.dump_pdb(f'BM3_{template.sequence()[mutation_position-1]}{mutation_position}{mutate_to}.pdb')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template',help= 'PDB template')
    parser.add_argument('-m', '--mutation',help= 'Format 82F')
    args = parser.parse_args()

    main(args)
