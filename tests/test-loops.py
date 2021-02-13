import  enz

def main():
    p = enz.protein('../data/1jme.pdb')
    enz.folds.init()
    pose = enz.folds.getpose(p.struc)
    print(enz.folds.fold_ccd(pose, {44:'Q'}))
    # enz.folds.fold_repack_mutate(pose, {44:'Q'})


if __name__ == '__main__':
    main()
