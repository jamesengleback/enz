# Enzyme-design
Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5750396/). ```enz``` is a thin wrapper for some structure prediction methods in ```pyrosetta``` and the autodock vina functionality in ```oddt```. ```enz``` is for enzyme design via rounds of mutant structure prediction  and substrate docking.

# Usage
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

# install
**From the installation folder**

```conda env create -f env.yml # create virtual environment with dependencies```

```conda activate enz```

```pip install . # install enz```


# todo
- loop & flexible region modelling
