import enz
import pandas as pd 
from tqdm import tqdm

from sequences import bm3_wt

def main():
    df = pd.read_csv('test_data/ManchesterSmiles.csv',
            index_col=0)

    p = enz.protein('../data/4key.pdb', 
            seq = bm3_wt,
            key_sites = [82, 87, 330],
            cofactors = ['HEM'])

    failures = []
    for i in tqdm(df['Smiles']):
        try:
            p.dock(i, exhaustiveness = 2)
        except Exception as e:
            failures.append(i)
    pd.Series(failures).to_csv('failures.csv')
                


if __name__ == '__main__':
    main()
