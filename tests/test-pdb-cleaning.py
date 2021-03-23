import enz 

def main():
    p = enz.protein('../benchmark/data/2UWH.pdb', cofactors=['HEM', 'PLM'])
    print(p.df['residue_name'].unique())
    p.save('test.pdb')

if __name__ == '__main__':
    main()
