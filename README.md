# Enzyme-design
Based on [**Small-molecule ligand docking into comparative models with Rosetta (2013)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5750396/), where the authors validate a method of structure prediction an ligand docking against known protein structures.  

# aim
- ```enz``` aims to provide a high level interface to some structure prediction methods from ```pyrosetta``` and docking methods from Autodock VINA

# 2 minute guide
```python
import enz

wt = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKE\
ACDESRFDKNLSQAWKFVRDFAGDGLVTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQ\
KWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKSQRANPDDP\
AYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIA\
GHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFS\
LYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRAC\
IGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK'

p = enz.protein('1jme.pdb',  # pdb path
                  seq = wt, # optional - aligns to structure
                  key_sites = [78,82,87,330, 181, 188]) # optional - constrain docking to this region

p.mutate(87,'V') # any residue
p.mutate(330, 'I')

p.refold() # repack side chains

p.save('new_structure.pdb')

results = p.dock('CCCCCCCCCCCC=O') # returns VINA score pandas.DataFrame

p.save_docking_results('docking_results') # save .pdb structures in new dir docking_results
```

# install from command line

```bash
git clone https://github.com/UoMMIB/enz.git
cd enz
conda env create -f env.yml
conda activate enz
pip install .
```
