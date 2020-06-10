# BM3-Design-PyRosetta
# I ❤️ BM3
# Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://github.com/jamesengleback/BM3-Design-PyRosetta/blob/master/docs/rosetta-ligand-dock-2013.pdf)

# Aim - BM3 structure prediction and Docking Tool
* **Renumbering** - the ```.pdb``` templates have missing residues, so they are aligned to the full sequence and map mutations in the sequence to the structure.
* **Folding**  - side chain packing is most important, also loop modelling for flexible regions (82-92, 435-439).
* **Docking** - autodock and rosetta docking are options. **Rosetta docking:** better performance **autodock:** 2-3 times faster.

 [**Computational methods and tools to predict cytochrome P450 metabolism for drug discovery**](https://www.ncbi.nlm.nih.gov/pubmed/30471192?otool=igbumllib) talks about using docking and other methods to predict P450 hydroxylation sites.
# progress
* ```/tools/tools.py - mutate_pose(pose, mm)``` mutates a pose of bm3 accoring to the mutations in ```mm``` -str; format: "a82f 87v l181K"  -not case sensitive, only mutate-to necessary. **packs side chains within 5A of mutation**
* ```quick_mutate.py``` - mutate bm3 template ```-t``` with mutations ```-x``` and save as ```-o``` usage:
```bash
python quick_mutate.py -t data/clean/3ben_clean.pdb -x "a82f f87v l181k i263p" -o output.pdb
```

# Folders
* **data** - BM3 fastas, ```.pdbs``` and clean ```.pdbs``` where water and ligands are removed (except heme).
* **testing** - testing folding methods
* **tools** - useful bits # todo - put renumerator here
* **tmp** - temporary
* **folding** - folding methods
* **docking** - autodock is easiest - try first
* **docs** - documentation for packages + papers

# Packages
I've saved an environment ```environment.yml``` run ```conda env create -f environment.yml``` to copy it, ```conda activate bm3``` to activate, ```conda deactivate``` to deactivate.

[**Pyrosetta**](http://www.pyrosetta.org/dow) - folding and docking
```
conda install pyrosetta -c  https://levinthal:paradox@conda.graylab.jhu.edu
```

[**scikit-bio**](http://scikit-bio.org/) - alignment
```
conda install -c https://conda.anaconda.org/biocore scikit-bio
```

[**ODDT - open drug discovery toolkit**](https://github.com/oddt/oddt)  paper [**docs**](https://oddt.readthedocs.io/en/latest/)
```
conda install -c oddt oddt
```

rdkit - dependency for oddt + useful cheminformatics tools
```
conda install rdkit -c rdkit
```

autodock-vina - docking
```
conda install -c bioconda autodock-vina
```
biopandas - cleaning pdbs
```
conda install -c conda-forge biopandas
```

# todo

## folding
The mutation programs I have re-pack sidechains with 5A of the mutation, additional options are:
- more sidechain repacking
- loop remodelling (82-92, 435-439)
- score functions - compared a bit in /data/stats
- validating models - backbone rmsd to targets?

## docking
- set up autodock script with odt - nn/rf score
- rosetta docking

# bits
- look into score outliers from /data/stats/stats.ipynb ; Group1 = 1smj, 3qi8, 3ekb, 3cbd ; Group2 = 3ekd, 4dqk, 4dql
