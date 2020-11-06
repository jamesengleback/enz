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
Install pyrosetta - download [here](http://www.pyrosetta.org/dow) - requires a "username" and "password", which are sent to you when yyou apply for a license. Make sure you get the right version for your machine. It shouldn't matter too much what version of python the pyrosetta distribution is for. Note that the download is really slow. On macOSX / linux, pyrosetta is distrubuted as a ```.tar.bz2``` archive. You can decompress these like this:
```bash
tar xfvj pypy3.7-v7.3.2-linux64.tar.bz2
```
then install by ```cd```'ing into ```PyRosetta4.Release.python37.linux.release-269/setup``` and running:
```bash
pip install .
```
```cd``` it where you want to download ```enz``` to
Clone enz to your machine:
```bash
git clone https://github.com/UoMMIB/enz.git
cd enz
```

