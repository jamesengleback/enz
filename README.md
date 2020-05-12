# BM3-Design-PyRosetta
# THE PLAN
**Throw science at the wall and see what sticks**

## Long Term
#### Aim
We want to predict roughly how different mutants of BM3 will interact with whatever substrate we give them. [**Computational methods and tools to predict cytochrome P450 metabolism for drug discovery.**](https://www.ncbi.nlm.nih.gov/pubmed/30471192?otool=igbumllib) mentions a couple of methods, including docking. For our job, we'll need:
* Structure Prediction
* Docking
* A good way to score/interpret the docking
* A program that pipes mutation|structure prediction|docking to find a useable mutant

## Short Term
### Structure Prediction
**Aim:** We want a program that (quickly 🤞 ) predicts the structure of a mutant.

I think we can rule out *ab-initio* folding for now because it takes so long! Plus we have plenty of structures to work from.

I'm not sure what the right method would be, but we could narrow it down to a handful and test all of them.

**Testing:** For now, I think a sensible way to test our algorithms is to try to predict known strucutures, and score based on that. alpha carbon (backbone) root mean squared devation (RMSD) has a function built into pyrosetta, but it won't take into account the side chains, so we'll need to think about that 🤔.

**Format of our scripts:** It would be good if the folding part of our scripts are contained in a function, so then we can import them into other scripts and test 'em.
