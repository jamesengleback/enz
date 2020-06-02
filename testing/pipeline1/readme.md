First attempt at a folding evaluation pipeline

### to do
* add folding methods to Folds.py, input format: fold(pose), no return object

### Aim:
* evaluate the accuracy and speed of different folding algorithms on known bm3 structures
* This pipeline will only start from BM3 WT, I'm considering starting from mutants in the future

### Steps
* randomly select a mutant bm3 structure ✅
* renumber the structures sequence and find mutations relative to wt
* make pose from wt
* copy wt pose and make mutations -> mutant pose
* fold - import any method from Folds.py
* score - ca rmsd + side chain positions between predicted and actual mutant structure


### Progress
* Mapping mutations onto template - just about works, one shortcoming is that it can't handle neighboring mutations, so has limited usability. In future, I might try a different alignment tool that doesn't do gap characters at mutations
* I've been randomly sampling from the clean pdb structures (using a relative path in the script) and renumbering according to the template. Should renumber according to wt sequence
* Folding - need to find a way to import folding methods, for now, just using backrub as a placeholder.
* Scoring - ca_rmsd is easy enough, need to find way to score sidechain similarity.

### Ideas
* Argparser - choose what folding method with shortcuts

### Questions
* **Renumbering:** should we renumber according to wt sequence or whatever's in the template structure - interprability thing.
alignments - see here http://scikit-bio.org/docs/0.2.0/alignment.html
