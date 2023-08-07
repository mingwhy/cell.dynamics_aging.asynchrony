"""
https://gitlab.com/olgaibanez/decibel/-/blob/main/module/decibel.py
Nov 17, 2022
"""  
def gcl(adata, num_divisions):
    """
    Following the original GCL.m script provided by the authors
    (https://github.com/guy531/gcl). 
    """
    from sklearn.metrics.pairwise import euclidean_distances
    import pandas as pd
    import numpy as np

    n, num_genes = adata.shape #changed
    #data = adata.X.todense().transpose() #changed
    data = adata.transpose() #changed by ming due to my input is a matrix
    gcl_output = np.zeros(num_divisions)

    def Vn(Aij, Bij):                    
        return 1/(n*(n-3)) * (np.sum(Aij*Bij) - n/(n-2)*np.sum(np.diag(Aij)*np.diag(Bij)))

    def Rn(Aij, Bij):
        return Vn(Aij, Bij)/np.sqrt(Vn(Aij, Aij)*Vn(Bij, Bij))
    
    def compute_Aij(X, n):
        d1 = euclidean_distances(X, X)
        m1 = np.mean(d1, 0)
        M1 = np.mean(d1)
        Aij = d1 - m1.T * np.ones((1, n))
        Aij = Aij.T - np.ones((n, 1)) * m1 
        Aij += M1
        Aij = Aij - d1 / n
        Aij[np.diag_indices(len(Aij))] = m1 - M1
        Aij = (n / (n-1)) * Aij
        return Aij


    for i in range(num_divisions):
        random_genes = np.arange(num_genes)
        np.random.shuffle(random_genes)

        X1 = data[random_genes[:num_genes//2],:].T
        X2 = data[random_genes[num_genes//2:],:].T
    
        Aij1 = compute_Aij(X1, n)
        Aij2 = compute_Aij(X2, n)
        
        gcl_output[i] = Rn(Aij2, Aij1)
        
    return gcl_output
