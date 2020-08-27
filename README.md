# Enzyme-design
Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5750396/), where the authors validate a method of structure prediction an ligand docking against known protein structures.  

### aim
```enz``` provides a simple interface to protein structure prediction and docking.

```enz``` works like this:
* structure prediction - ```pyrosetta``` - side chain repacking within 5 A of mutation site ✅ flexible region remodelling (cyclic coordinate descent) - in progress ⏳
* ligand docking - autodock VINA via ```oddt``` ✅
* scoring of docking poses using VINA ✅ KD prediction using nnscore - in progress ⏳



### Description
🔨 under construction! 🔨 

### ```enz.protein```
 ```enz.protein```: contains the protein structure and interfaces with ```pyrosetta```. ```enz.protein``` can be initialised using a protein sequence as well as a structure, which is aligned to structure, allowing for canonical numbering or mutation be inserting new FASTA sequences. ```enz.protein``` automaticly cleans ```pdb``` structures by removing duplicated chains, water and non-cofactor ligands.

 ```python
enz.protein(pdb_path, seq=None) # initial params, seq optional
```

 ```enz.protein``` has 3 important methods:
 * ```.mutate(<aa_num (int)>, <aa (str, capital)>)``` which mutates the sequence, ready for refolding
 * ```.refold``` which finds all diffences between the desired sequence and the structure's sequence, and mutates the residues with side chain repacking within 5 A of the mutation site
 * ```.save(<save file path>)``` saves the predicted structure as a ```.pdb``` file

### ```enz.vina```
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

### simple mutate & save

```python
import enz

bm3WT= 'MTIKEMPQPKTFGELKNLP...'


p = enz.protein(pdb_path = '1jme.pdb', seq = bm3WT) # optional: align sequence for canonical numbering

p.mutate(82, 'F')
p.mutate(87, 'V')
p.refold()
p.save('A82F-F87V.pdb')
```

### screen compound library

```python
import enz
import pandas as pd
import os

df = pd.read_csv('fda.csv') # fda compound library

bm3WT= 'MTIKEMPQPKTFGELKNL...'

p = enz.protein('3ben.pdb', seq=bm3WT) # align a sequence for canonical numbering

vina = enz.vina(p, active_site_aas = [400, 82, 87, 263, 188, 400, 330], exhaustiveness=16)
# initialise a vina object, narrow search space to some active site amino acids


for name, smiles in zip(df['Drug Name'], df['SMILES']):
    results= vina.dock(smiles,name)
    scores = results.autodock_score # results also contains coordinates of poses
    if os.path.exists('fda-screen-results.csv'):
	results.autodock_score.to_csv('fda-screen-results.csv', mode='a', header=False) # add to file (append mode)
    else:
	results.autodock_score.to_csv('fda-screen-results.csv') # new file (first run)
```
### alanine scan & dock with one substrate


# installation guide (terminal)

* clone enz: ```$ git clone https://github.com/UoMMIB/enz.git```
* move into enz: ```$ cd enz```
* create  virtual environment with dependencies: ```$ conda env create -f env.yml```
* activate: ```$ conda activate enz```
* install enz: ```$ pip install . ```


# todo
- loop & flexible region modelling
- cache cleaning
- nn score
