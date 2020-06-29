from enz import  enz, tools
import pandas as pd
import os
from tqdm import tqdm

def clean(path):
    # clean pdb
    cache = '__cleaning-cache__'
    os.makedirs(cache, exist_ok=True)
    pdb_data = tools.clean_pdb(path)
    out_filename = os.path.join(cache, os.path.basename(path))
    pdb_data.to_pdb(out_filename)
    return out_filename

def main():
    seq = tools.fasta_to_series('data/sequences/A8-4-prenyltransferase-aa.fasta')[0]
    struc = clean('data/raw/lewis/5kcg.pdb')
    pt = enz.Protein(pdb_path = struc, seq = seq)
    pt.mutate_seq(256, 'L')
    pt.mutate_seq(399, 'A')
    pt.mutate_seq(17,'T')
    pt.refold()
    pt.dump('pt-mutant.pdb')

if __name__ == '__main__':
    main()
