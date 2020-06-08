# BM3-Design-PyRosetta
# I ❤️ BM3
# Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://github.com/jamesengleback/BM3-Design-PyRosetta/blob/master/rosetta-ligand-dock-2013.pdf)

# Aim - BM3 structure prediction and Docking
## Structure Prediction
* **Renumbering** - the ```.pdb``` templates have missing residues, so they are aligned to the full sequence and map mutations in the sequence to the structure.
* **Folding**  - side chain packing is most important, also loop modelling for flexible regions (82-92, 435-439).
* **Docking** - autodock and rosetta docking are options. **Rosetta docking:** better performance **autodock:** 2-3 times faster.

 [**Computational methods and tools to predict cytochrome P450 metabolism for drug discovery**](https://www.ncbi.nlm.nih.gov/pubmed/30471192?otool=igbumllib) talks about using docking and other methods to predict P450 metabolites.

 # Folders
* **data** - BM3 fastas, ```.pdbs``` and clean ```.pdbs``` where water and ligands are removed (except heme).
* **testing** - testing folding methods
* **tools** - useful bits # todo - put renumerator here
* **tmp** - temporary
* **folding** - folding methods
* **docking** - docking
