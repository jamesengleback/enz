#!/bin/python
import re
import pandas as pd


def extract(text):
    text = text.split('\n')
    table_start = ['---+--' in i for i in text].index(True) + 1
    table = []
    for row in text[table_start:]:
        items = row.split()
        is_all_ints = lambda l : sum([re.search('-?\d+', i) is not None for i in l])
        if len(items) == 4 and is_all_ints(items):
            table.append(dict(zip(['mode','affinity (kcal/mol)', 'dist from best mode - rmsd - ub','dist from best mode - lb'], items)))
    print(pd.DataFrame(table))


def main():
    with open('output.txt','r') as f:
        data = f.read()
    print(extract(data))

if __name__ == '__main__':
    main()
