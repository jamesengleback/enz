import pandas as pd
import skbio
import re
from tqdm import tqdm
import os
import random
import argparse
import pyrosetta
import Folds


def FastaToSeries(path):
    # opens fasta files, outputs series
    # full fasta string in one cell
    with open(path,'r') as file:
        data = file.read()
    data = [i.split('\n') for i in data.split('>')]
    index = [i[0] for i in data][1:] # first item is ''
    values = [''.join(i[1:]) for i in data][1:] # first item is []
    df = pd.Series(values, index=index)
    return df

def Renumber(query_seq, return_both=False):
    # canonical numbering
    ref_seq = '''MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRL\
    IKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
    KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYD\
    ENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTS\
    GLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVL\
    GGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATL\
    VLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRKKGC'''.replace(' ','')
    # sequences need to be in the skbio.Protein object
    ref_seq = skbio.Protein(ref_seq)
    query_seq = skbio.Protein(query_seq) # item unpacks the string
    # align - slow, gives efficiency warning - skbio devs working on this
    # check for updates in future
    print('\x1b[31m') #red
    aln = skbio.alignment.global_pairwise_align_protein(ref_seq,query_seq)
    print('\x1b[0m') # reset colors
    # unpack aligned sequences
    aln_ref = ''.join([i.decode() for i in aln[0].loc[0].values])
    aln_query = ''.join([i.decode() for i in aln[0].loc[1].values])
    if return_both:
        return dict(enumerate(aln_query,0)), dict(enumerate(aln_ref,0))
    else:
        return dict(enumerate(aln_query,0))


def MapMutations(query_seq, mutations):
    # 1. renumber query_seq - map to original
    # 2. index mutations - use same alignment
    # 3. return mutation positions in
    d_query, d_ref = Renumber(query_seq, return_both=True)
    tmp = {i:d_query[i] for i in d_query if d_query[i] != '-'}
    mapping = {j:i for i,j in enumerate(tmp, 1)} # original:renumbered, start count at 1
    mutations = {}
    for i, j in zip(d_query, d_ref):
        query_aa = d_query[i]
        ref_aa = d_ref[j]
        if query_aa != ref_aa and query_aa != '-' and ref_aa != '-' and ref_aa != 'Z':
            mutations[i] = {'from':ref_aa, 'to':query_aa}
    mapped_mutations = {mapping[i]:mutations[i]['to'] for i in mutations}
    return mapped_mutations

def ParseMutation(m):
    number = int(''.join([i for i in m if i.isdigit()]))
    letter = [i for i in m if i.isalpha()][-1] # in case they give 2 letters
    return number, letter.upper()

def MakeMutations(pose,mutations):
    for i in mutations:
        print('\x1b[31m') #red
        print(i,mutations[i])
        pyrosetta.toolbox.mutate_residue(pose, i,mutations[i])

def main(args):
    #pyrosetta.init()
    # get template
    template_file = '../../data/clean/1jme_clean.pdb' # wt???
    print('\x1b[31m') #red
    print(f'Template = {os.path.basename(template_file)}')
    #template = pyrosetta.pose_from_pdb(template_file)

    mutations = dict([ParseMutation(i) for i in args.mutations.split()])
    print(mutations)

    # make mutant from template
    #mutant_pose = pyrosetta.Pose()
    #mutant_pose.assign(template)

    #mutations = MapMutations(mutant_pose.sequence())
    #MakeMutations(mutant_pose, mutations)
    #for i in range(0,len(template.sequence())-50,50):
    #    print()
    #    print(template.sequence()[i:i+50])
    #    print(mutant_pose.sequence()[i:i+50])


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
