# Enzyme-design
# Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://github.com/jamesengleback/BM3-Design-PyRosetta/blob/master/docs/rosetta-ligand-dock-2013.pdf)

# Aim - Homology structure prediction and Docking Tool
# Usage
```python
import enz

enzyme = enz.protein(pdb = '****.pdb', seq = 'mtikemlmssp...')

enzyme.mutate(pos = 82,aa= 'F') # lazy
# or
enzyme.mutate(seq = seq) # lazy

enzyme.fold(loop = True) # generate ensemble

vina = enz.vina(protein = enzyme,
                cpds = 'c1ccccc', # or iterable
                cpd_names = 'benzene', # same len() as cpds
                score = None,
                output_dir = None)

def custom_score(protein, cpd):
  ...
  return score     

vina.score = custom_score

vina.dock() # dock ensemble, return score

```
