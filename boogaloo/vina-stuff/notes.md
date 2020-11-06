# enz upgrade
## aims
- remove bloat - enz relies on oddt for its autodock vina capability only, and on skbio for its sequence alignment.
- vina - use subprocess
- alignment - implement needleman wunsch
- pyrosetta - don't rely on anaconda channel - use local install.


# vina.py
- vina class - inits with ligand smiles and protein pdb
  - why class and not function?
- makes tempfile with pdbqt structures
- uses 10 as all center + sizes
  - todo: residue active site selection
- parse output


# vina-fn.py
- vina docking as a function
- clean pdb
  - removes water
  - chooses a chain (A by default)
  - remove ligands (except for cofactors - e.g. heme)
