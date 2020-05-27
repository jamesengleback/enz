First attempt at a folding evaluation pipeline

### Aim:
* evaluate the accuracy and speed of different folding algorithms on known bm3 structures
* This pipeline will only start from BM3 WT, I'm considering starting from mutants in the future

### Steps
* randomly select a mutant bm3 structure
* renumber the structures sequence and find mutations relative to wt
* make pose from wt
* copy wt pose and make mutations -> mutant pose
* fold - any method
* score - ca rmsd + side chain positions between predicted and actual mutant structure


### Progress
* Mapping mutations onto template - just about works, one shortcoming is that it can't handle neighboring mutations, so has limited usability. In future, I might try a different alignment tool that doesn't do gap characters at mutations
* Folding - need to find a way to import folding methods, for now, just using backrub as a placeholder.
* Scoring - ca_rmsd is easy enough, need to find way to score sidechain similarity.
