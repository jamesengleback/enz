from biopandas.pdb import PandasPdb
import os
import sys
import argparse
from tqdm import tqdm

colors = {
 'Reset': '\x1b[0m',
 'BrightBlack': '\x1b[30;1m',
 'BrightRed': '\x1b[31;1m',
 'BrightGreen': '\x1b[32;1m',
 'BrightYellow': '\x1b[33;1m',
 'BrightBlue': '\x1b[34;1m',
 'BrightMagenta': '\x1b[35;1m',
 'BrightCyan': '\x1b[36;1m',
 'BrightWhite': '\x1b[37;1m'}

def GetSeq(file):
    data = PandasPdb().read_pdb(file)
    seq = data.amino3to1()['residue_name'].to_list()
    seq = ''.join(seq)
    return seq

def main(args):
    '''
    If input is a folder, try to get sequence from all of them
    else Just the one
    If output is none, print sequence
    '''
    input = args.input

    sequences = {}
    if input != None:
        if os.path.isdir(input):
            # if folder
            files = os.listdir(input)
            # check for .pdb extension
            files = [i for i in files if '.pdb' in i]
            print(colors['BrightGreen'])
            for i in tqdm(files):
                path_to_file = os.path.join(input, i)
                sequences[f'>{i}'] = GetSeq(path_to_file)

        #if single file
        else:
            sequences[os.path.basename(input)] = GetSeq(input)
    else:
        print(colors['BrightRed'], 'No input')
        print('Exiting')
        sys.exit


    output = args.output
    if output != None:
        if os.path.isdir(output):
            output = os.path.join(output, 'Sequences.fasta')
        with open(output,'w') as f:
            for i in sequences:
                f.write(i + '\n')
                # split into lines of 60
                seq=sequences[i]
                splitline_sequence = [seq[j:j+60] + '\n' for j in range(0,len(seq), 60)]
                f.writelines(splitline_sequence)
                f.write('\n')
    else:
        print(colors['BrightBlue'])
        for i in sequences:
            print(i)
            print(sequences[i])



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='input file or folder')
    parser.add_argument('-o','--output',help='output path')
    args = parser.parse_args()
    main(args)
