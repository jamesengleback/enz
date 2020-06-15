# tools
# ```tools.py```
* ```fasta_to_series(path)``` -read fasta file as pandas series
* ```aln(s1,s2)``` - align s1 and s1, return: strings of aligned s1 and s2
*  ``` diff(s1,s2)``` - align, find mutation in s2 relative to s1, return: dict: ```{1:{'from':'A','to':'P'}}```
* ```mutate_sequence(seq, pos,aa)``` - simple mutate position in string - 1 indexed, return:mutated string
* ```map_sequences(s1,s2)``` - align sequences, return: dict of {num_s1:num_s2} mapping of sequence index
* ```parse_mutation(m)``` - m: str format "a82f"/"82f" (not case sensitive);  returns: (82, "F")
* ```parse_mutations(mm)``` - mm: str format eg. "a82f f87v l181k ", slipts and returns parse_mutation(m) for m in mm
* ```mutate_pose(pose, mm)``` - pose & mm as above, mapping handled internally, repack radius = 5.0 A returns None, acts inplace

# demo
I think ```mutate_pose``` is most useful, needs to be tested properly, inside, I set the sidechain repack radius to 5A.
```python
import tools
import pyrosetta

template_file = '../data/clean/3ben_clean.pdb'
mm = "394H 395H 396H 397H 398H 399H 401H 402H 403H" # around C400
outfile = 'mutant.pdb'
pyrosetta.init()
bm3 = pyrosetta.pose_from_pdb(template_file)
tools.mutate_pose(bm3, mm) # make mutants
bm3.dump_pdb(outfile)
```
