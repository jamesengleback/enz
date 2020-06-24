# Enzyme-design
# Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://github.com/jamesengleback/BM3-Design-PyRosetta/blob/master/docs/rosetta-ligand-dock-2013.pdf)

# Aim - Homology structure prediction and Docking Tool
# Usage
```python
import enz

path = 'data/clean/1jme_clean.pdb'
wt = tools.fasta_to_series('bm3-wt.fasta')[0]
bm3 = enz.Protein(pdb_path = 'data/clean/1jme_clean.pdb', seq = wt)

bm3.mutate_seq(82, 'F')
bm3.mutate_seq(87, 'V')
bm3.refold()

# todo: bm3.dock(...)
vina = enz.Vina(bm3)
scores, results = vina.dock('c1ccccc', 'benzene')
# scores = pd.DataFrame
# results = [oddt.mol, ...]
```

# install

```conda env create -f env.yml```

```conda activate enz```

```python setup.py install```


# todo
- setup.py
- oddt rdkit instead of obabel
- test on more examples
- env
- nnscore
