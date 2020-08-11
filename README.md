# Enzyme-design
Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5750396/), where the authors validate a method of structure prediction an ligand docking against known protein structures.  

### aim
```enz``` provides a simple interface to protein structure prediction and docking.

```enz``` works like this:
* structure prediction - ```pyrosetta``` - side chain repacking within 5 A of mutation site ✅ flexible region remodelling (cyclic coordinate descent) - in progress ⏳
* ligand docking - autodock VINA via ```oddt``` ✅
* scoring of docking poses using VINA ✅ KD prediction using nnscore - in progress ⏳



### Description
 ```enz``` has two object types:

### ```enz.protein```
 ```enz.protein```: contains the protein structure and interfaces with ```pyrosetta```. ```enz.protein``` can be initialised using a protein sequence as well as a structure, which is aligned to structure, allowing for canonical numbering or mutation be inserting new FASTA sequences. ```enz.protein``` automaticly cleans ```pdb``` structures by removing duplicated chains, water and non-cofactor ligands.

 ```python
enz.protein(pdb_path, seq=None) # initial params, seq optional
```

 ```enz.protein``` has 3 important methods:
 * ```.mutate_seq(<aa_num (int)>, <aa (str, capital)>)``` which mutates the sequence, ready for refolding
 * ```.refold``` which finds all diffences between the desired sequence and the structure's sequence, and mutates the residues with side chain repacking within 5 A of the mutation site
 * ```.dump(<save file path>)``` saves the predicted structure as a ```.pdb``` file

### ```enz.Vina```
```enz.Vina``` wraps the ```oddt``` vina interface, which launches a VINA docking simulation. ```enz.Vina``` can be initialised from a ```.pdb``` file or from an ```enz.protein``` object.

```python
# initial params
enz.Vina(
    protein=None, # enz.protein or .pdb file (clean?)
    center=None, # 3D coords of active site (optional)
    box_dims=None,# tuple xyz size of simulation box (optional)
    acitve_site_aas=None, # auto-draw box around active site !!!!!
    ncpus=None, # default = ncpus -1
    exhaustiveness=8,# number of poses to predict
)
```

Behind the scenes, VINA sets up a box within the strucutre to simulate within, but it helps to narrow it down to your active site. If you know your active site amino acid numbers, then initialise like: ```enz.Vina(<enz.protein>, acitve_site_aas = [1,2,3,4, ... ])```

```enz.Vina``` has one main methods:
* ```enz.dock(<smiles str>, <cpd name (optional)>)``` which starts a docking run. Currently, it returns a ```pandas.DataFrame``` with the Autodock scores (ref?) for each pose, and a list of ```oddt.mol``` objects.

# examples

### simple mutate & dock benzene
```python
import enz

p = enz.protein(pdb_path = '1jme.pdb') # optional: align sequence

p.mutate_seq(82, 'F')
p.mutate_seq(87, 'V')
p.refold()
p.dump('A82F-F87V.pdb') # save

vina = enz.Vina(p) # init vina from enz.protein

scores, results = vina.dock('c1ccccc', 'benzene') # scores: pd.DataFrame; results: [oddt.mol, ...] (poses)
scores.to_clipboard()
```
### screen compound library

```python
df = pd.read_csv('fda.csv') # fda compound library

bm3WT= 'MTIKEMPQPKTFGELKNLPLLNTDKPVQA\
LMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFV\
RDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMV\
RALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGE\
QSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALY\
FLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTA\
PAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPER\
FENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTN\
YELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK*'

p = enz.protein('3ben.pdb', seq=bm3WT) # align a sequence for canonical numbering

vina = enz.Vina(p, acitve_site_aas = [400, 82, 87, 263, 188, 400, 330], exhaustiveness=16) # initialise a vina object, narrow search space to some active site amino acids

# save an empty dataframe with same format as 'scores'
pd.DataFrame([],columns = ['vina_affinity', 'vina_rmsd_lb', 'vina_rmsd_ub', 'vina_rmsd_input',
   'vina_rmsd_input_min', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion',
   'vina_hydrophobic', 'vina_hydrogen', 'SMILES Atom Order', 'cpd']).to_csv('fda-screen-results.csv')

for name, smiles in zip(df['Drug Name'], df['SMILES']):
    try:
        scores, poses = vina.dock(smiles,name)
        scores.to_csv('fda-screen-results.csv', mode='a', header=False) # add to file (append mode)
    except:
        print('🤦') # handle errors
```
### alanine scan & dock with one substrate


# installation guide (terminal)
### clone repository:
Make sure you have git installed
```$ git clone https://github.com/UoMMIB/enz.git```
### move into ```enz``` folder:
```$ cd enz```

### create & activate virtual environment with dependencies (takes some time, make sure you have good internet connection):
```$ conda env create -f env.yml```
### activate:
```$ conda activate enz```
### install enz:
```$ pip install . ```


# todo
- loop & flexible region modelling
