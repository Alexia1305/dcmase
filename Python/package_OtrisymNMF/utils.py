import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.preprocessing import normalize
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import chi2_kernel


def bestMap(L1, L2):
    L1 = np.array(L1).flatten()
    L2 = np.array(L2).flatten()
    if L1.shape != L2.shape:
        raise ValueError("L1 and L2 must have the same size")
    Label1, Label2 = np.unique(L1), np.unique(L2)
    L1_new = np.array([np.where(Label1 == l)[0][0] for l in L1])
    L2_new = np.array([np.where(Label2 == l)[0][0] for l in L2])
    nClass = max(len(Label1), len(Label2))
    G = np.zeros((nClass, nClass), dtype=int)
    for i in range(len(L1)):
        G[L1_new[i], L2_new[i]] += 1
    row_ind, col_ind = linear_sum_assignment(G.max() - G)
    mapping = {col: row for row, col in zip(row_ind, col_ind)}
    newL2_new = np.array([mapping[l] for l in L2_new])
    newL2 = np.array([Label1[l] for l in newL2_new])
    return newL2

def accuracy(L_true, L_pred):
    L_map = bestMap(L_true, L_pred)
    return 100.0 * np.sum(L_true == L_map) / len(L_true)

def build_multigraph(Person, sigma=1.0):
    X = Person.T
    dist_euc = cdist(X, X, metric='euclidean')

    # noyaux
    X_rbf = np.exp(-dist_euc ** 2 / (2 * sigma ** 2))
    X_cauchy = 1.0 / (1.0 + (dist_euc ** 2 / sigma ** 2))
    X_chi2 = chi2_kernel(X, gamma=1.0 / (2 * sigma ** 2))
    X_inter = np.minimum(X[:, None, :], X[None, :, :]).sum(axis=2)

    X_multi = np.stack([X_rbf, X_cauchy, X_chi2, X_inter], axis=0)
    return X_multi