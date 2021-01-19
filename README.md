# Enzyme-design
enz wraps protein structure prediction methods from pyrosetta and molecular docking from autodock vina for template-based enzyme design. 

# features
## alignment
```enz``` can align an amino acid sequence to a structure to handle your preferred numbering system. any differences between the input sequence and the structure sequences will be mutated. 
## structure cleaning
```enz``` automatically cleans ```.pdb``` structures by removing water, duplicated chains and other small molecules that aren't specified as cofactors.
## folding
enz uses ```pyrosetta``` for structure prediction functions. currently ```enz``` only repacks side-chains around the mutation site. planned feature: loop remodelling for flexible regions.
## docking
```enz``` uses ```autodock vina``` for molecular docking, since it's fairly fast and straightforward. file conversion is handled automatically and the target area is defined by numbering residues. results are scored using vina's built-in scoring system. results objects, like protein objects, contain ```pandas``` ```DataFrames``` of atomic coordinates which can be used to generate custom score functions. results objects can also be saved fairly easily.  currently, side chains are treated as rigid. a planned feature is to enable side chain flexibility in the target site. a known issue is that some atom types (e.g. boron) are rejected by vina. I'll see what i can  do about that.


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

results = p.save('docking_results') # save .pdb structures in new dir docking_results
```

# install from command line
## conda
You'll need ```conda```. I've made an environment file that you can automatically create the  environment from with the command
```python
conda create -f env.yml # execute from the enz directory```

Install ```pyrosetta``` - download [here](http://www.pyrosetta.org/dow) - requires a "username" and "password", which are sent to you when yyou apply for a license. Make sure you get the right version for your machine. It shouldn't matter too much what version of python the ```pyrosetta``` distribution is for. Note that the download is really slow. On macOSX / linux, ```pyrosetta``` is distrubuted as a ```.tar.bz2``` archive. You can decompress these like this:
```bash
tar xfvj pypy3.7-v7.3.2-linux64.tar.bz2
```
then install by ```cd```'ing into ```PyRosetta4.Release.python37.linux.release-269/setup``` and running:
```bash
pip install .
```
## other requirements
### openbabel
### nwalign3
### biopandas

```cd``` to where you want to download ```enz``` to
Clone enz to your machine:
```bash
git clone https://github.com/UoMMIB/enz.git
cd enz
```

