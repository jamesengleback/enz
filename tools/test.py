'''
test tools.py
'''
import tools
import pyrosetta


def main():
    # args
    template_file = '../data/clean/3ben_clean.pdb'
    #mm = "a82f f87v l181k i263p"
    #mm = "11H 12H 13H 14H 15H 16H 17H 18H 19H a82H 83H 84H 85H 86H 87H"
    mm = "394H 395H 396H 397H 398H 399H 401H 402H 403H" # around C400
    outfile = 'mutant.pdb'
    pyrosetta.init()
    bm3 = pyrosetta.pose_from_pdb(template_file)
    tools.mutate_pose(bm3, mm)
    bm3.dump_pdb(outfile)

if __name__ == '__main__':
    main()
