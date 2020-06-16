def find_heme(protein):
    atoms = [i.atomicnum for i in list(protein.atoms)]
    coords = protein.coords
    d = {idx:{'atom':a, 'coords':c} for idx, a, c in zip(range(len(atoms)), atoms, coords)}
    heme_fe = [d[i] for i in d if d[i]['atom'] == 26]
    return next(iter(heme_fe)) # gets it out of a list

def euclidian_distance(a,b):
    return np.linalg.norm(a-b)

def distance_to(mol,target):
    d = {}
    for atom in mol:
        a = np.array(atom.coords)
        d[atom.idx0] = euclidian_distance(a, target['coords'])
        # extract more info
    print(d)
