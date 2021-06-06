# benchmark
## aim
Assess accuracy of docking in ```enz``` by mutating a template structure to X and docking ligand, docking accuracy.
## data
```data/``` contains a collection of ligand-bound P450 BM3 mutants. ```strucs.csv``` contains the corresponding ligand SMILES, its name in the pdb file and the mutations for each structure
## template selection
base paper uses blast to find close structure, I might just ro random perms.
## score
- rmsd cutoff
