# todo
- dev: run hrac experiment
- dev: try new enzymes
- Task: Package
- share: intro session: plan
- Task: loop modelling - sahara
- Task: docking score - NN-Score
- Task: docking score - custom
- Task: validate docking - pdb
- Task: validate structure
- Task: multiprocess
- Task: detect active site + narrow box
- plan: MIB-github - talk to Sam + linus


# autodock score error 20200629
Run - screen hrac with 1 structure
```
Parse error on line 14 in file "/tmp/autodock_vina_2rtf0408/ligands_3n95xvsf/0.pdbqt": ATOM syntax incorrect: "As" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.
 33%|████████████████████▋                                         | 98/293 [1:54:03<3:46:57, 69.83s/it]
Traceback (most recent call last):
  File "hrac-screen.py", line 33, in <module>
    main()
  File "hrac-screen.py", line 29, in main
    scores, results = vina.dock(s,n)
  File "/home/jamese/miniconda3/envs/enz/lib/python3.8/site-packages/enz/enz.py", line 140, in dock
    scores = self.autodock_score(results, self.vina)
  File "/home/jamese/miniconda3/envs/enz/lib/python3.8/site-packages/enz/enz.py", line 122, in autodock_score
    scores = pd.concat([pd.Series(i.data) for i in vina.score(results)],
  File "/home/jamese/miniconda3/envs/enz/lib/python3.8/site-packages/pandas/core/reshape/concat.py", line 271, in concat
    op = _Concatenator(
  File "/home/jamese/miniconda3/envs/enz/lib/python3.8/site-packages/pandas/core/reshape/concat.py", line 329, in __init__
    raise ValueError("No objects to concatenate")
ValueError: No objects to concatenate
```

- Offender = ```DSMA    C[As](=O)([O-])[O-]```
- Error on As ```ATOM syntax incorrect: "As" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.```
- similar: http://autodock.1369657.n2.nabble.com/ADL-Atom-types-td4028847.html " add to parameters in the parameter file for ligand"
- solution: **add atom types to parameter file**
- add new atom types in autodock: http://autodock.scripps.edu/faqs-help/faq/how-do-i-add-new-atom-types-to-autodock-4
- quickfix = ```try: .. except: ..```
