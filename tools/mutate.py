import pandas as pd
import skbio
import re
from tqdm import tqdm
import os
import random
import argparse
import pyrosetta
import Folds

def ParseMutation(m):
    if m != None:
        number = int(''.join([i for i in m if i.isdigit()]))
        letter = [i for i in m if i.isalpha()][-1] # in case they give 2 letters
        return number, letter.upper()
    else:
        return None

def MutateSequence(mutations):
    # mutations is list
    # return mutated sequence
    # map to pose later
    seq = list('''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ',''))
    for num, aa in mutations:
        seq[num] = aa
    return ''.join(seq)


def MutatePose(mutant, pose):
    print('\x1b[31m') #red
    mutant_seq = skbio.Protein(mutant)
    pose_seq = skbio.Protein(pose.sequence())
    aln = skbio.alignment.global_pairwise_align_protein(mutant_seq,pose_seq)
    aln_mutant = ''.join([i.decode() for i in aln[0].loc[0].values])
    aln_pose = ''.join([i.decode() for i in aln[0].loc[1].values])
    new_pose_seq = []
    for i,j in zip(aln_mutant, aln_pose):
        if i != j and j != '-' and i != '-' :
            new_pose_seq.append(i)
        else:
            new_pose_seq.append(j)
    new_pose_seq = ''.join(new_pose_seq)
    new_pose_seq = new_pose_seq.replace('-','')
    pose_seq_ascii = ''.join([i.decode() for i in pose_seq.values] )
    mutations = [{'index':idx,'from':i,'to':j} for idx,(i,j) in enumerate(zip(pose_seq_ascii, new_pose_seq),1) if i != j]
    for i in mutations:
        idx, old, new = i['index'], i['from'], i['to']
        pyrosetta.toolbox.mutate_residue(pose, idx,new)
        print(f' \x1b[31m mutation: {old + str(idx) + new}')
    print('\x1b[0m') # reset colors



def test():
    #new_seq = MutateSequence([ParseMutation(i) for i in args.mutations.split()])
    template_file = args.template
    #all_templates = [os.path.join('../data/clean/', i) for i in os.listdir('../data/clean/')]
    #template_file = random.sample(all_templates, 1)[0]

    pyrosetta.init()
    template_pose = pyrosetta.pose_from_pdb(template_file)
    # make mutant from template
    mutant_pose = pyrosetta.Pose()
    mutant_pose.assign(template_pose)
    MutatePose(new_seq, mutant_pose)

    print(f'template: \t {os.path.basename(template_file)}')
    diff(mutant_pose.sequence(), template_pose.sequence()) #check, takes time though

    # fold
    Folds.PackChainsMini(mutant_pose)

    if args.output == None:
        mutant_pose.dump_pdb(f"../tmp/BM3-{args.mutations.replace(' ','-')}.pdb")
    else:
        mutant_pose.dump_pdb(args.output)

if __name__ == '__main__':
    test()
