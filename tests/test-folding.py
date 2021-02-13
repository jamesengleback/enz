import enz

def main():
    p = enz.protein('../data/4key.pdb')

    p.mutate(44, 'Q')

    p.refold()

    print(p.df)

if __name__ == '__main__':
    main()
