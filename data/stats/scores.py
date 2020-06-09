import pandas as pd
from tqdm import tqdm
import os
import pyrosetta

def score(path):
    name = os.path.basename(path).split('_')[0]
    pose = pyrosetta.pose_from_pdb(path)
    fa = pyrosetta.pyrosetta.get_fa_scorefxn()
    # PyRosetta4.Release.python37.ubuntu.release-253/pyrosetta/database/scoring/weights
    st = pyrosetta.create_score_function('score0') # weights tag - database weights file
    d= {'full atom':fa(pose)}
    #{'standard':st(pose)}
    return pd.Series(d, name = name)

def main():
    pyrosetta.init()
    clean_pdbs = [os.path.join('../clean/',i) for i in os.listdir('../clean') if 'pdb' in i]
    df = pd.DataFrame([score(i) for i in clean_pdbs])

    print(df)
if __name__ == '__main__':
    main()
