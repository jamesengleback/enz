## Useful bits

### Renumbering
In Mutations.ipynb and renumber.py I almost made a tool that renumbers our sequences to canonical numbering. It uses alignments, because the structure sequences are missing segments.

The issues I was having was that the alignment tool would put in a gap at every mismatch, so I've used regex to find letter-letter patterns and remove the spacer on both the template and query alignments. This is just a quick fix though, and wouldn't work if there were two neighboring mutations. It doesn't work for 1SMI.pdb.
