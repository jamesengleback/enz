from openbabel import pybel 

def main():
    smiles = 'CCCC=O'
    m = pybel.readstring('smi',smiles)
    m.addh()
    m.make3D()
    m.write('pdbqt', 'test.pdbqt', overwrite=True)

if __name__ == '__main__':
    main()
