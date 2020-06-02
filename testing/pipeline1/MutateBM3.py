import pandas as pd
from Bio import pairwise2
import re
from tqdm import tqdm
import os
import random
import argparse

import pyrosetta

import Folds


def renumber(query_seq, return_both=False):
    seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLS\
    SQRLIKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYH\
    AMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVR\
    ALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKD\
    PETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSY\
    KQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWG\
    DDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYEL\
    DIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')

    alns = pairwise2.align.globalxx(seq,query_seq,penalize_end_gaps = False,
    one_alignment_only=True)
    s1, s2, _,_,_ = alns[0]
    s1 = re.sub('\w-\w','',s1)
    s2 = re.sub('\w-\w','',s2)
    if return_both:
        return dict(enumerate(s1,0)), dict(enumerate(s2,0))
    else:
        return dict(enumerate(s2,0))

def ParseMutation(m):
    number = int(''.join([i for i in m if i.isdigit()]))
    letter = [i for i in m if i.isalpha()][-1] # in case they give 2 letters
    return number, letter.upper()

def MakeMutations(pose,mutations):
    # forward renumber
    pose_sequence = pose.sequence()
    d1 = renumber(pose_sequence)
    d1 = {i:d1[i] for i in d1 if d1[i] != '-'}
    print(d1)
    for (num, aa) in mutations:
        d1.update({num:aa})
    s1 = ''.join([d1[i] for i in d1])
    for i in range(5):
        ix = i*100
        print(pose_sequence[ix:ix+100])
        print(s1[ix:ix+100])
        print()




def main(args):
    pyrosetta.init()

    # get template
    template_file = '../../data/clean/1jme_clean.pdb' # wt???
    print('\x1b[31m') #red
    print(f'Template = {os.path.basename(template_file)}')
    template = pyrosetta.pose_from_pdb(template_file)

    # make mutant from template
    mutant_pose = pyrosetta.Pose()
    mutant_pose.assign(template)

    d1,d2 = renumber(mutant_pose.sequence(), return_both=True)

    print('\x1b[31m') #red
    mutations = [ParseMutation(i) for i in args.mutations.split(' ')]
    MakeMutations(mutant_pose, mutations) #fix


    # imported them, access with this dictionary
    #FoldMethods = {'backrub':Folds.BackrubMinimise,
    #'sidechains':Folds.PackSideChainsSimple}
    #fold = FoldMethods[args.fold]
    #print('\x1b[31m') #red
    #fold(mutant_pose)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutations')
    parser.add_argument('-f', '--fold')
    args = parser.parse_args()
    main(args)
