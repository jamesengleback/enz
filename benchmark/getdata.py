import os
import io
import requests 
import pandas as pd
from tqdm import tqdm

def main():
    if not os.path.exists('data'):
        os.makedirs('data')
    df = pd.read_csv('strucs.tsv', delimiter='\t')
    for i, j in tqdm(zip(df['code'], df['link']), total = len(df)):
        try:
            r = requests.get(j)
            if r.status_code == 200:
                data = r.text
                with open(os.path.join('data', f'{i.upper()}.pdb'), 'w') as f:
                    f.write(data)
        except Exception as e:
            print(e)


if __name__ == '__main__':
    main()
