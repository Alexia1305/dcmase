import numpy as np
from sklearn.preprocessing import normalize
from scipy.linalg import orthogonal_procrustes, svd, orth
from scipy.sparse import csr_matrix
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
def clustering_QR_embeddings(W):
    L, n, r = W.shape

    # Binariser la matrice W
    for l in range(L):
        Wl = W[l]
        W_bin = np.zeros_like(Wl)
        W_bin[np.arange(n), np.argmax(Wl, axis=1)] = 1.0
        W[l] = W_bin

    #
    Q_ref = W[0].copy()
    for l in range(1, L):
        R, _ = orthogonal_procrustes(W[l], Q_ref)
        W[l] = W[l] @ R

    coherence = np.zeros(L)
    for i in range(L):
        for j in range(L):
            if i != j:
                coherence[i] += np.linalg.norm(W[i].T @ W[j], 'fro')
    weights = coherence / np.sum(coherence)

    Z = np.hstack([
        np.sqrt(weights[l]) * W[l]
        for l in range(L)
    ])

    U, S, _ = np.linalg.svd(Z, full_matrices=False)
    k_opt = min(r, np.searchsorted(np.cumsum(S ** 2) / np.sum(S ** 2), 0.98) + 1)
    Q = U[:, :k_opt]

    Q -= np.mean(Q, axis=0, keepdims=True)
    Q_norm = normalize(Q, norm='l2', axis=1)

    return Q_norm

def nmf_consensus(labels, r):
    L, n = labels.shape

    rows, cols, data = [], [], []
    col = 0
    for l in range(L):
        for c in range(r):
            idx = np.where(labels[l] == c)[0]
            rows.extend(idx)
            cols.extend([col] * len(idx))
            data.extend([1] * len(idx))
            col += 1

    H = csr_matrix((data, (rows, cols)), shape=(n, L*r))
    Q = NMF(n_components=r, init="nndsvd").fit_transform(H)
    return Q.argmax(axis=1)





def consensus_kmeans(labels, k, normalize_rows=True):
    """
    Consensus clustering via binaire matrix + KMeans.

    Parameters
    ----------
    labels : ndarray of shape (L, n)
        Cluster assignments of each point in each partition.
    k : int
        Number of consensus clusters.

    normalize_rows : bool
        Normalize rows of H to unit L2 norm .

    Returns
    -------
    consensus_labels : ndarray of shape (n,)
        Cluster assignment for each point.
    """


    L, n = labels.shape


    # Construction de la matrice binaire H (dense)
    H = np.zeros((n, int(L * k)), dtype=float)
    col_offset = 0
    for l in range(L):
        for c in range(k):
            idx = np.where(labels[l] == c)[0]
            H[idx, col_offset] = 1.0
            col_offset += 1

    # Normalisation des lignes (optionnelle)
    if normalize_rows:
        row_norms = np.linalg.norm(H, axis=1, keepdims=True) + 1e-12
        H = H / row_norms

    # KMeans dense
    kmeans = KMeans(n_clusters=k, random_state=42)
    consensus_labels = kmeans.fit_predict(H)

    return consensus_labels

