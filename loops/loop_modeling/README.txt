loop modeling script

This script does loop modeling over residues 77-85 on the PDB file
"test_in.pdb".  Loop modeling uses a combination small/shear moves, fragment
insertion, and CCD loop closure, in a simulated annealing MC protocol to
generate closed, low-energy loop conformations. Packaged with this script are
the input PDB file ("test_in.pdb") and the corresponding 3-residue fragment
file ("test_in.frag3")

Instructions:

1) to generate a fragment library for your PDB of choice, go to
http://robetta.bakerlab.org/fragmentsubmit.jsp and enter the sequence of the
PDB file in FASTA format.
2) once you have the fragment file, edit the script to use your desired input
files and define the beginning and end residue #'s of your loop
3) run the script.

Note: you must create a fragment file for each new PDB protein sequence that
you use.
 
