from biopandas.pdb import PandasPdb
import os
import argparse
import sys
from tqdm import tqdm

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
 'BrightWhite': '\x1b[37;1m',
 'BackgroundBlack': '\x1b[40m',
 'BackgroundRed': '\x1b[41m',
 'BackgroundGreen': '\x1b[42m',
 'BackgroundYellow': '\x1b[43m',
 'BackgroundBlue': '\x1b[44m',
 'BackgroundMagenta': '\x1b[45m',
 'BackgroundCyan': '\x1b[46m',
 'BackgroundWhite': '\x1b[47m',
 'BackgroundBrightBlack': '\x1b[40;1m',
 'BackgroundBrightRed': '\x1b[41;1m',
 'BackgroundBrightGreen': '\x1b[42;1m',
 'BackgroundBrightYellow': '\x1b[43;1m',
 'BackgroundBrightBlue': '\x1b[44;1m',
 'BackgroundBrightMagenta': '\x1b[45;1m',
 'BackgroundBrightCyan': '\x1b[46;1m',
 'BackgroundBrightWhite': '\x1b[47;1m',
 'Bold': '\x1b[1m',
 'Underline': '\x1b[4m',
 'Reversed': '\x1b[7m'}

def clean_pdb(pdb_path):
    data = PandasPdb().read_pdb(pdb_path)
    singleChainData = SingleOutChain(data)
    strippedDown = StripHetAtoms(singleChainData)
    return singleChainData


def SingleOutChain(pandas_pdb):
    # return complete pandas pdb object
    unique_chains = pandas_pdb.df['ATOM']['chain_id'].unique()
    chain_1 = unique_chains[0]
    pandas_pdb.df['ATOM'] = pandas_pdb.df['ATOM'].loc[pandas_pdb.df['ATOM']['chain_id'] == chain_1]
    return pandas_pdb

def StripHetAtoms(pandas_pdb):
    # strips water and ligands, except the heme
    pandas_pdb.df['HETATM'] = pandas_pdb.df['HETATM'].loc[pandas_pdb.df['HETATM']['residue_name'] == 'HEM']
    return pandas_pdb


def main(args):
    # file stuff
    input = args.input
    output = args.output
    if input == None:
        print(colors['BrightRed'], 'No input specified!')
        print('Specify an input PDB file or folder with -i path/to/input')
        sys.exit()
    input_parent_folder =  os.path.dirname(os.path.dirname(input))
    default_output_folder = os.path.join(input_parent_folder, 'clean')

    if output == None:
        print(colors['BrightYellow'], 'No output file specified')
        if os.path.exists(default_output_folder):
            # default folder already exists
            # make new one
            print(f'Making default folder: {default_output_folder}')
            print('Coud be overwriting exisitng folder ;)')
            os.makedirs(default_output_folder, exist_ok = True)
            output = default_output_folder

    if os.path.isdir(input):
        files = os.listdir(input)
        files = [i for i in files if '.pdb' in i]
        for i in tqdm(files):
            pdb_code = i.split('.')[0]
            biopandas_obj = clean_pdb(os.path.join(input,i))
            out_file_name = f'{pdb_code}_clean.pdb'
            biopandas_obj.to_pdb(os.path.join(output, out_file_name))

    if os.path.isfile(input):
        # need to make compatible with single files
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help = 'input file or folder')
    parser.add_argument('-o', '--output', help = 'output file or folder')
    args = parser.parse_args()
    main(args)
