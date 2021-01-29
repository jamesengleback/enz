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
## dataframes
all molecule objects have a ```mol.df``` proterty, which displays a pandas dataframe of the molecule's ```.pdb```, including its coordinates which can be very useful for creating custom docking scores.


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

p.save('new_structure.pdb') # save any moelcule object as a pdb with the save() method

results = p.dock('CCCCCCCCCCCC=O') # returns a results object
# which contains the docking poses as enz.mol objects and 
# the calculated binding energy of each

results = p.save('docking_results') # save .pdb structures in new dir docking_results
```

# install from command line
## 1. clone enz
if you have ```git``` then clone this repository to your machine
```sh
git clone https://github.com/UoMMIB/enz.git
```
get in the repo!
```sh
cd enz
```

## 2. conda environment setup
You'll need ```conda```. I've made an environment file that you can automatically install most of the dependencies. set it up with:

```sh
conda create -f env.yml # execute from the enz directory
```

then activate it with 

```sh
conda activate enz
```
you'll need to activate this environment before using enz. 

# 3. install pyrosetta in the ```enz``` environment
Install ```pyrosetta``` - download [here](http://www.pyrosetta.org/dow) - requires a "username" and "password", which are sent to you when you apply for a license. Make sure you get the right version for your machine's operating system. I've set up the environment for python 3.7, so best get the python 3.7 version. Note that the download is really slow. On macOSX / linux, ```pyrosetta``` is distrubuted as a ```.tar.bz2``` archive. You can decompress these like this:
```bash
tar xfvj pypy3.7-v7.3.2-linux64.tar.bz2
```
then install by ```cd```'ing into ```PyRosetta4.Release.python37.linux.release-269/setup``` and running:
```bash
pip install .
```

## 4. back to enz
Navigate back to ```enz``` and install with
```sh
pip install .
```
at the base of the ```enz``` file tree

that's it. ```enz``` is a pain in the ass to install because of ```pyrosetta```, but for now there's not much i can do about that.

# how to use enz for enzyme design:
## ```enz.protein```
You'll need a ```pdb``` file template of your protein to work on. Optionally, the amino acid sequence and the names of the cofactors as the occur in the ```pdb``` file. Initialising an ```enz.protein``` object automatically cleans the ```pdb``` structure by removing duplicated chains, water and any other molecules not specified in the ```cofactors = [...``` argument. Cleaning the protein in this way is necessary for rosetta and vina compatibility. 
You might want to include an amino acid sequence if you'd like to use a residue numbering system that differes from that of the ```pdb``` file. The sequence you provide will be aligned to that of the structure and any differences between it and the sequence of the structure will be resolved at ```refold()``` 

```python
import enz

sequence = 'MSAKBNGFUIAUIEA...'

p = enz.proten('XXXX.pdb', # essential
		cofactors = ['NADP'], # optional, must be a list
		seq = sequence) # optional
```
## dataframes
```enz.protein``` and ```enz.mol``` objects (docking results) both have a ```df```  property, which gives a ```pandas.DataFrame``` of the coordinates of each atom. This includes the ```x_coord``` ```y_coord``` & ```z_coord``` for each atom, which is useful to know if you want to score docking on the basis of how close two atoms are together, for example. 

```python
p.df

>>>  record_name  atom_number blank_1  ... element_symbol charge line_idx
0           ATOM            1          ...              N    NaN      645
1           ATOM            2          ...              C    NaN      646
2           ATOM            3          ...              C    NaN      647
3           ATOM            4          ...              O    NaN      648
4           ATOM            5          ...              C    NaN      649
...          ...          ...     ...  ...            ...    ...      ...
```

## ```mutate()```  & ```refold()```
Mutate positions in the protein with the ```mutate(<position int>, <amino acid letter str>)```. The amino acid letters are case insensitive. Nothing is actually calculated until you call the ```refold()``` method. ```refold()``` replaces sidechains in the structure at all positions where the aligned ```enz.protein.seq``` differs from that of the structure. side chains are repacked within a radiums of 5 Angstroms of the mutation site by default, but this can be tweaked with ```refold(10)``` for example. This method calls ```pyrosetta``` which spits out a huge ammount of text.

```python
p.mutate(55, 'A')

p.mutate(99,'P')
for i in [100, 122, 188]:
	p.mutate(i, 'A')

p.refold()
```
## ```save()```
the ```enz.protein.save('filename.pdb')``` saves a ```.pdb``` file of the protein to a desired location. 
## ```dock()```
dock ligands to your protein with the ```enz.protein.dock(...``` method. The method requires the SMILES code for your compound and a list of ```target_residues``` - a list of amino acid positions (int) around which the simulation box will be drawn. If you want to have the compounds bind to the active site, then provide some numbers of residues that are in the active site. you don't need all of them, just enough to draw a box around.

an optional argument to ```dock(``` is ```exhaustiveness``` - use an integer 1-16, for low-high resolution docking respectively.

the ```dock()``` method returns a results object that wraps the poses and the calculated affinity of each, as well as a ```save()``` method.

```python
results = p.dock('CCCCCCCCCC=O', target_residues = [50, 55, 80, 199, 330])
```
## results - ```enz.vina.results```
contains the docking results poses and a ```DataFrame``` of the scores. 
### attributes:
- ```.scores``` - ```pandas.DataFrame``` of each ligand's binding energy as calculated by vina. ```mode``` refers to the ligand id.
- ```dictionary``` provides access to each docking pose and their binding energy individually.
#### save
- the ```enz.vina.results``` object has a ```save(...) ``` method saves the docking results into a new directory specified as an argument. the directory contains a ```csv``` of the scores dataframe, the receptor as a ```pdb``` and a ```pdb``` of each pose

