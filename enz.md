# enz

# improve
-  enz.Protein(path, seq=None)    
- active site detection?? light weight
- docstrings!!!!!!!!!
- nnscore - lower priority
- rdkit backend????

# tests
# master, 20200728:
# In [14]: bm3 = enz.Protein(path, seq=s)                                                                                                                                  
# In [15]: bm3.mutate_seq(82,'F')                                                                                                                                          
# In [16]: bm3.refold()                                                                                                                                                    
- behaves: seq aln AFTER refold ???
# In [19]: bm3.dump('~/ops/')                                                                                                                                              
- core.io.pdb.pdb_writer: {0} [ ERROR ] StructFileRep::dump_pdb: Unable to open file:~/ops/bm3.pdb for writing!!!
- FileNotFoundError: [Errno 2] No such file or directory: '~/ops/bm3.pdb'
- works in current directory


# munro
- benchling protocols e.g. tev  - get from alessia, type up; share pdfl; make #protocols
