import enz
import numpy as np
import pandas as pd
from tqdm import tqdm # progress bar

BM3_WT = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK*'

def score(protein, results):
    fe = protein.df.loc[protein.df['atom_name']=='FE',['x_coord', 'y_coord', 'z_coord']]
    distances, energies = [], []
    for i,j in zip(results.df['mode'], results.df['affinity (kcal/mol)']):
        p_i = results.poses[i] # enz.mol object
        c1 = p_i.df.loc['atom_number' == 1,['x_coord', 'y_coord', 'z_coord']]
        distances.append(np.linalg.norm(fe - c1))
        energies.append(j)
    d, e = np.array(distances), np.array(energies)
    # scale e 0:1
    e -= e.min()
    e /= e.max()
    # mean weighted distance by binding energy
    x = e * d
    return x.mean()


def main():
    df = pd.DataFrame([], columns=['mutation','score'])

    for i in tqdm([75,87,263,181,188]):
        p = enz.protein('../data/4key.pdb', seq = BM3_WT, cofactors = ['HEM'])
        r = p.dock('CCCCCC=CCC=CCC=CCC=CCCCC(=O)O', target_residues=[400,49,181,330,262]) # arachidonic acid, active site
        r.save(f'bm3_{p.seq[i]}{i}A') # eg bm3_L75A
        df.append(pd.DataFrame([i, score(p, r)], columns = ['mutation','score']))
    df.to_csv('alanine-scan.csv', index=False)
if __name__ == '__main__':
    main()

