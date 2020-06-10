# BM3-Design-PyRosetta
# I ❤️ BM3
# Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://github.com/jamesengleback/BM3-Design-PyRosetta/blob/master/docs/rosetta-ligand-dock-2013.pdf)

# Aim - BM3 structure prediction and Docking Tool
* **Renumbering** - the ```.pdb``` templates have missing residues, so they are aligned to the full sequence and map mutations in the sequence to the structure.
* **Folding**  - side chain packing is most important, also loop modelling for flexible regions (82-92, 435-439).
* **Docking** - autodock and rosetta docking are options. **Rosetta docking:** better performance **autodock:** 2-3 times faster.

 [**Computational methods and tools to predict cytochrome P450 metabolism for drug discovery**](https://www.ncbi.nlm.nih.gov/pubmed/30471192?otool=igbumllib) talks about using docking and other methods to predict P450 hydroxylation sites.

# Folders
* **data** - BM3 fastas, ```.pdbs``` and clean ```.pdbs``` where water and ligands are removed (except heme).
* **testing** - testing folding methods
* **tools** - useful bits # todo - put renumerator here
* **tmp** - temporary
* **folding** - folding methods
* **docking** - autodock is easiest - try first
- **docs** - documentation for packages + papers

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
OpenBabel - mol file handling
```
conda install -c openbabel openbabel
```

# todo

## folding
- sidechain repacking - restrict to active site
- loop remodelling (82-92, 435-439)
- score functions - which

## docking
- set up autodock script with odt
- rosetta docking

## bits
- save conda env + leave instructions
- move tesing tools into tools
- clean testing
- link https papers
- compare score functions
