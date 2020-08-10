# Enzyme-design
Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5750396/). ```enz``` is for enzyme design via rounds of mutant structure prediction  and substrate docking.

### Description
```enz``` is a thin wapper for simple structure prediction methods in ```pyrosetta``` and substrate docking methods in ```autodock vina``` facilitated by ```oddt```. ```enz``` has two object types:

### ```enz.protein```
 ```enz.protein```: contains the protein structure and interfaces with ```pyrosetta```. ```enz.protein``` can be initialised using a protein sequence as well as a structure, which is aligned to structure, allowing for canonical numbering or mutation be inserting new FASTA sequences. ```enz.protein``` automaticly cleans ```pdb``` structures by removing duplicated chains, water and non-cofactor ligands.

 ```python
enz.protein(pdb_path, seq=None) # initial params, seq optional
```

 ```enz.protein``` has 3 mportant methods:
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

### alanine scan & dock with one substrate


### screen compound library




# install
In the terminal:
###
### create & activate virtual environment with dependencies
```conda env create -f env.yml```

```conda activate enz```
### install enz
```pip install . ```


# todo
- loop & flexible region modelling
