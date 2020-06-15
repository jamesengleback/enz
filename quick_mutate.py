import sys
sys.path.append('tools')
import tools
import pyrosetta
import argparse


def main(args):
    # input params
    TEMPLATE = args.template
    MM = args.mutations
    OUT = args.output

    # control flow + default params
    if TEMPLATE==None:
        print('No template PDB (-t /path)')
        sys.exit()
    if MM == None:
        MM = ""
    if OUT == None:
        OUT='out.pdb'

    # protein stuff
    pyrosetta.init()
    bm3 = pyrosetta.pose_from_pdb(TEMPLATE)
    tools.mutate_pose(bm3, MM)
    bm3.dump_pdb(OUT)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template')
    parser.add_argument('-x', '--mutations')
    parser.add_argument('-o','--output')
    args = parser.parse_args()
    main(args)
