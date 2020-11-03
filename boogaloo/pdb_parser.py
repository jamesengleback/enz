
class PDB:
    def __init__(self, path):
        self.path = path
        self.data = self.read()
    @property
    def df(self):
        pass
    def read(self):
        with open(self.path, 'r') as f:
            data = f.readlines()
        return data
    


def main():
    pdb = PDB('test_data/4key.pdb')
    print(pdb.data)
if __name__ == '__main__':
    main()
