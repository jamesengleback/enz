import pandas as pd
from tqdm import tqdm
import os
import pyrosetta

def score(path):
    name = os.path.basename(path).split('_')[0]
    pose = pyrosetta.pose_from_pdb(path)
    fa = pyrosetta.pyrosetta.get_fa_scorefxn()
    # PyRosetta4.Release.python37.ubuntu.release-253/pyrosetta/database/scoring/weights
    s0 = pyrosetta.create_score_function('score0') # weights tag - database weights file
    s1 = pyrosetta.create_score_function('score1')
    s2 = pyrosetta.create_score_function('score2')
    s3 = pyrosetta.create_score_function('score3')
    s5   = pyrosetta.create_score_function('score5')
    ref = pyrosetta.create_score_function('ref2015')
    dk = pyrosetta.create_score_function('docking')
    d= {'full atom':fa(pose),
    #'score0':s0(pose),
    #'score1':s1(pose),
    #'score2':s2(pose),
    #'score3':s3(pose),
    #'score5':s5(pose),
    'ref2015':ref(pose),
    'docking':dk(pose)
    }
    return pd.Series(d, name = name)

def main():
    pyrosetta.init()
    clean_pdbs = [os.path.join('../clean/',i) for i in os.listdir('../clean') if 'pdb' in i]
    raw_pdbs = [os.path.join('../raw/',i) for i in os.listdir('../raw') if 'pdb' in i]
    df_clean = pd.DataFrame([score(i) for i in clean_pdbs])
    df_raw = pd.DataFrame([score(i) for i in raw_pdbs])
    print(df_clean)
    print(df_raw)
    df_raw.to_csv('BM3-rosetta-scores-raw.csv')
    df_clean.to_csv('BM3-rosetta-scores-clean.csv')
if __name__ == '__main__':
    main()
