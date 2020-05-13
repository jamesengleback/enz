### data/sequences
I used ```tools/GetSequences.py ``` to extract the sequences, and then muscle with default parameters:
```muscle -in Sequences.fasta -out Sequences_msa.fasta```
Installing muscle with homebrew, try one of:
```bash
brew install homebrew/science/muscle
brew install homebrew/science/netcdf
brew install homebrew/science/dssp
```
Linux:
```bash
sudo apt-get install muscle
```

### To Do:
 * Filter out structures/sequences that aren't BM3 Heme domain
 * Find the mutations of each sequence


 #### Potential solutions:
 ```python
 In [1]: import pandas as pd                                                                                               

 In [2]: df = pd.read_csv('Sequences_msa.fasta', lineterminator='>',header=None)                                           

 In [3]: df = df[0].str.split('\n',expand=True)                                                                            

 In [4]: df_2 = df[1] # make a fresh copy of first column                                                                                                     

 In [5]: df_2.index = df[0] # index with pdb names
 ### to do: join all the columns into a single string
 ### then split each character into its own cell
 ### then drup duplicates to find each position
 ```
