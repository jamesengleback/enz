import numpy as np

score = lambda a, b : 1 if a == b else 0

def score_matrix(s1,s2):
    mat = np.zeros(shape = (len(s1) + 1, len(s2) + 1), dtype=int)
    mat[0] = range(0,-len(mat[0]),-1)
    mat[:,0] = range(0,-len(mat[:,0]),-1)
    # weird indexing here, room for errors
    for rownum, row in enumerate(mat[1:],1):
        for colnum, col in enumerate(row[1:],1):
            s1_i, s2_j = s1[rownum-1], s2[colnum-1]
            match_score = score(s1_i, s2_j)
            indel_score = -1
            left, above, diag = mat[rownum, colnum-1] + indel_score,\
                    mat[rownum-1,colnum] + indel_score,\
                    mat[rownum-1,colnum-1] + match_score
            cell_score = max(left, above, diag)
            mat[rownum, colnum] = cell_score
    return mat

def traceback(score_mat, s1, s2):
    pass

def main():
    s1 = 'CGATCGATGCCCG'
    s2 = 'CGTAGCTAGGTCTA'
    print(score_matrix(s1,s2))
    
if __name__ == '__main__':
    main()
