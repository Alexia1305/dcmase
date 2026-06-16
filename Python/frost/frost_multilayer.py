
import time
from scipy.sparse import issparse, find, csr_matrix
from scipy.sparse.linalg import norm
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import normalized_mutual_info_score

from frost.aggregation_way.co_association import clustering_coassociation, MatrixFreeSpectralClusteringCoAssociation, \
    USENC_ConsensusFunction
from frost.aggregation_way.QR_embeddings import clustering_QR_embeddings, nmf_consensus, consensus_kmeans
import numpy as np
import math
from .SVCA import svca
from scipy.sparse import diags


# import Cluster_Ensembles as CE


# -----------------------------------------------
# ---------------- frost ----------------
# -----------------------------------------------

def frost_multilayer(X_list, r, numTrials=3, maxiter=1000, delta=1e-6, time_limit=300, init_method='MF-SC-CA', init_w=None,
                     power_method=False, init_partition=None, verbosity=0,
                     init_seed=None, true_labels=None):
    """
    Heuristic algorithm for multilayer community detection via joint nonnegative matrix trifactorization.
    Estimates nonnegative matrices S_l>=0 and W_l>=0 that minimize:

        sum_{l=1..L} ||X_l-W_l S_l W_l^T||_F^T

    subject to the constraints:
        W_l = D_l V
        W_l^T W_l = I
        S_l >= 0
        W_l >= 0
    where:
    - X_l is the adjacency matrix of layer l
    - S_l is a nonnegative community interaction matrix for layer l
    - W_l is a nonnegative structured orthogonal matrix encoding community memberships
    - D_l is a diagonal scaling matrix specific to layer l
    - V is a shared binary community assignment matrix across all layers

    Notes
    -----
    The matrices W_l=D_l V are not stored explicitly but represented with:

    - v : ndarray of shape (n,) — layer-independent community assignments defining V,
      where v[i] is the community index of node i.

    - w : ndarray of shape (L, n) — layer-dependent diagonal values defining D_l,
      where w[l, i] is the scaling factor for node i in layer l.

    Thus, each row i of W_l has a single nonzero element:

        W_l[i, v[i]] = w[l, i]

    The vector v is shared across layers, while w is layer-dependent.

     Parameters:
         X_list : list of crs_matrix
             List of L matrices, one per layer.
             - Each matrix has shape (n, n), where n is the number of nodes.
             - Matrices are nonnegative and symmetric for adjacency matrix of an undirected graphs.
         r : int
             Number of communities.
         numTrials : int, default=1
             Number of trials with different initializations.
         maxiter : int, default=1000
             Maximum iterations for each trial.
         delta : float, default=1e-7
             Convergence tolerance (Stop if error<delta or error_prec-error<delta).
         time_limit : int, default=300
             Time limit in seconds.
         init_method : str, default=SVCA
             Initialization method ("random", "SVCA").
         verbosity : int, default=1
             (1 for messages, 0 for silent mode).
         init_seed : float, optional (default=None)
             Random seed for the initialization for the experiments

     Returns:
         w_best : ndarray of shape (L, n)
            Layer-specific nonzero values of W_l (diagonal elements of D_l).
         v_best : np.array, shape (n,)
             Shared community assignment vector. v_best[i] is the community of node i.
         S_best : ndarray of shape (L, r, r)
             Layer-specific community interaction matrices S_l.
         error_best : float
             Relative error ||X - ZSZ'||_F / ||X||_F.
         time_per_iteration : list[list[float]]
             Time per iteration for each trial.
             Each inner list contains the times (in seconds) for every iteration
             of that specific trial.
     """
    errors = []
    start_time = time.time()
    # Only one layer must be in a list
    if (isinstance(X_list, np.ndarray) or issparse(X_list)) and X_list.ndim == 2 and X_list.shape[0] == X_list.shape[1]:
        X_list = [X_list]

    if not isinstance(X_list, list):
        raise TypeError("Xlist must be a Python list")

    X_list = [csr_matrix(X).astype('float') if not issparse(X) else X.tocsr().astype('float') for X in X_list]
    L = len(X_list)
    n = X_list[0].shape[0]
    error_best = float('inf')


    w_best = np.zeros((L, n))
    S_best = np.zeros((L, r, r))
    v_best = np.zeros(n)

    # Precomputations
    normX = [norm(X, 'fro') for X in X_list]
   
    degrees_layers = []

    for X in X_list:
        degrees = np.asarray(X.sum(axis=1)).ravel()
        degrees_layers.append(degrees)

    if verbosity > 0:
        print(f'Running {numTrials} Trials in Series')

    for trial in range(numTrials):
        if init_partition is not None:
            v=init_partition.copy()
            w = np.zeros((L, n))
            for l in range(L):
                w[l] = initialize_w_values(X_list[l], v)

        else:

            if L == 1 :
                init_method = 'onelayer'

            if init_seed is not None:
                init_seed += 10 * trial
                np.random.seed(init_seed)

            if init_method == 'random':
                base = np.arange(r)
                rest = np.random.randint(0, r, size=n - r)
                v = np.concatenate([base, rest])
                np.random.shuffle(v)
                w = np.zeros((L, n))
                for l in range(L):
                    w[l] = initialize_w_values(X_list[l], v)
                # for l in range(L):
                #     w[l], v = PowMethOTRISYMNMFFixed(X_list[l], r, labels_final, w[l], maxiter=50, timelimit=200)

            elif init_method == 'onelayer':
                w = np.zeros((L, n))
                lr = np.random.randint(0, L)  # Choose a random layer
                w[lr], v, _ = initialize_W_onelayer(X_list[lr], r)
                # default value degree of the node / sum degrees same community
                for l in range(L):
                    if l == lr:
                        continue
                    w[l] = initialize_w_values(X_list[l], v)
                # for l in range(L):
                #     w[l], v = PowMethOTRISYMNMFFixed(X_list[l], r, labels_final, w[l], maxiter=50, timelimit=200)


            else:

                w, v = initialize_W_alllayers(X_list, r, init_method, init_w, power_method)

        # Normalization of W
        for l in range(L):
            nw = np.zeros(r)
            for i in range(n):
                nw[v[i]] += w[l, i] ** 2
            nw = np.sqrt(nw)

            denom = nw[v]
            mask = denom != 0
            w[l, mask] /= denom[mask]
            w[l, ~mask] = 0

        # Compute of S
        S = update_S(X_list, r, w, v)

        if verbosity:
            if true_labels is not None:
                print('NMI initialisation : ', normalized_mutual_info_score(true_labels, v))
            print('Time', time.time() - start_time)
        prev_error = 0
        for l in range(L):
            prev_error += compute_error(normX[l], S[l])
        prev_error=prev_error/sum(normX)
        error = prev_error
        init_error = (error)

        errors.append(init_error)
        for iteration in range(maxiter):
            if time.time() - start_time > time_limit:
                print('Time limit passed')
                break

            w, v = update_W(X_list,degrees_layers, S, w, v)

            S = update_S(X_list, r, w, v)

            prev_error = error
            error = 0
            for l in range(L):
                error += compute_error(normX[l], S[l])
            error=error/sum(normX)
            

            if error < delta or abs(prev_error - error) < delta:
                break

        if error < error_best:
            w_best, v_best, S_best, error_best = (
                 w.copy(), v.copy(), S.copy(), error
            )

        if verbosity > 0:
            print(f'Trial {trial + 1}/{numTrials} with {init_method}: Error {error:.4e} | Best: {error_best:.4e}')
            if true_labels is not None:
                print('NMI : ', normalized_mutual_info_score(true_labels, v_best))
            print('Time', time.time() - start_time)

    return w_best, v_best, S_best, errors

def update_W(X_list, degrees_layers, S, w, v):
    L = len(X_list)
    n = X_list[0].shape[0]
    r = S.shape[1]

    """
    # Pre-calculations to avoid a double loop on ‘n’
    """
    wp2 = np.zeros((L, r))
    S2 = S**2
    w2 = w**2
    Xii = np.array([X.diagonal() for X in X_list])

    for l in range(L):
        for k in range(r):
            wp2[l, k] = np.sum(w2[l]*S2[l, v, k])

    # Update of each row (node)
    for i in range(n):
        vi_new = -1
        wi_new = np.full(L, -1)
        wi = np.zeros((L, n))
        f_new = np.inf


        # Test each community
        for k in range(r):
            wi = w.copy()
            vi = v.copy()
            vi[i] = k
            erreur = 0
            # For each layer, find the best value for w[i] with v[i] = k
            for l in range(L):
                if degrees_layers[l][i]==0 :
                    wi[l][i] = 0
                else:
                    c3 = S2[l, k, k]
                    c1 = 2 * (wp2[l, k] - (w[l, i] * S[l, v[i], k]) ** 2) - 2 * S[l, k, k] * Xii[l, i]
                    X = X_list[l]
                    start = X.indptr[i]
                    end = X.indptr[i + 1]

                    cols = X.indices[start:end]
                    vals = X.data[start:end]

                    mask = cols != i
                    selected_cols = cols[mask]
                    selected_vals = vals[mask]

                    c0 = -4 * np.sum(selected_vals * w[l, selected_cols] * S[l, v[selected_cols], k])

                    # Résolution des racines avec la méthode de Cardan
                    roots = cardan_depressed(4 * c3, 2 * c1, c0)

                    # Trouver la meilleure solution positive pour w_l(i, k_i)
                    x = 0
                    min_value = c3 * (x ** 4) + c1 * (x ** 2) + c0 * x
                    for sol in roots:
                        value = c3 * (sol ** 4) + c1 * (sol ** 2) + c0 * sol
                        if sol > 0 and value < min_value:
                            x, min_value = sol, value

                    wi[l][i] = x

                    erreur += min_value

            if erreur < f_new:
                f_new = erreur
                wi_new = wi[:, i].copy()
                vi_new = k

        for l in range(L):
            for k in range(r):
                wp2[l, k] = wp2[l, k] - (w[l, i] * S[l, v[i], k]) ** 2 + (wi_new[l] * S[l, int(vi_new), k]) ** 2

        # Mise à jour des valeurs de v et w
        v[i] = vi_new
        for l in range(L):
            w[l, i] = wi_new[l]

    # Normalization of W
    for l in range(L):
        nw = np.zeros(r)
        for i in range(n):
            nw[v[i]] += w[l, i] ** 2
        nw = np.sqrt(nw)

        denom = nw[v]
        mask = denom != 0
        w[l, mask] /= denom[mask]
        w[l, ~mask] = 0

    return w, v


def cardan_depressed(a, c, d, tol=1e-12):
    """ Cardano formula to find the roots of ax^3+cx+d=0 """
    if abs(a) < tol:
        if abs(c) < tol:
            return []
        else:
            return [-d / c]

    # b=0 t^3+pt+q
    p = c / a
    q = d / a
    Delta = 4 * (p ** 3) + 27 * (q ** 2)

    if abs(Delta) < tol:
        return [0]
    elif Delta > 0:  # one real solution
        sqrtD = np.sqrt(Delta / 27)
        return [np.cbrt((-q + sqrtD) / 2) + np.cbrt((-q - sqrtD) / 2)]

    else:  # 3 real different solutions or multiple solution

        r = 2 * np.sqrt(-p / 3)
        cos_arg = -q / 2 * np.sqrt(-27 / (p ** 3))
        cos_arg = np.clip(cos_arg, -1, 1)
        theta = np.arccos(cos_arg) / 3
        return r * np.cos(np.array([theta, theta + 2 * np.pi / 3, theta + 4 * np.pi / 3]))


def update_S(X_list, r, w, v):
    L = len(X_list)
    S = np.zeros((L, r, r))
    for l in range(L):

        # Get the row indices, column indices, and values of the non-zero elements in the sparse matrix X
        i, j, val = find(X_list[l])

        # Loop through the non-zero elements of X
        for k in range(len(val)):
            S[l, v[i[k]], v[j[k]]] += w[l, i[k]] * w[l, j[k]] * val[k]
    return S



# ------------------------------------------------
# -------------- INITIALIZATION ------------------
# ------------------------------------------------

def initialize_W_alllayers(X_list, r, init_method, init_w, power_method):
    L = len(X_list)
    n = X_list[0].shape[0]

    W = np.zeros((L, n, r))
    w = np.zeros((L, n))
    v_init = np.zeros((L, n))
    # Find w and v for each layer
    for l in range(L):
        w[l], v_init[l], W[l] = initialize_W_onelayer(X_list[l], r)

    if init_method == 'co-association':
        C = clustering_coassociation(W)  # calcul de la matrice de co-association à partir de W
        model = SpectralClustering(n_clusters=r, affinity='precomputed', random_state=42)
        labels_final = model.fit_predict(C)

    elif init_method == 'QR':
        C = clustering_QR_embeddings(W)
        clustering_model = KMeans(n_clusters=r, n_init=50, random_state=42)
        labels_final = clustering_model.fit_predict(C)

    elif init_method == 'NMF':
        labels_final = nmf_consensus(v_init, r)
    elif init_method == 'simple':
        labels_final = consensus_kmeans(v_init, r, normalize_rows=True)
    elif init_method == 'MF-SC-CA':
        method = MatrixFreeSpectralClusteringCoAssociation(v_init, r)
        labels_final = method.fit_predict()

    elif init_method == 'USENC':
        labels_final = USENC_ConsensusFunction(v_init.T, r)

    for l in range(L):
        w[l] = initialize_w_values(X_list[l], labels_final)

    # for l in range(L):
    #     w[l], v = PowMethOTRISYMNMFFixed(X_list[l], r, labels_final, w[l], maxiter=50, timelimit=200)

    return w, labels_final

def initialize_W_onelayer(X, r):
    options = {'average': 1}
    n = X.shape[0]
    p = max(2, math.floor(0.1 * n / r))
    # Estimation of ZS=Z*S with SVCA
    # Attention the rank of X must be greater than r
    ZS, K = svca(X.tocsc(), r, p, options=options)
    norm2x_squared = X.multiply(X).sum(axis=0)  # matrice 1 x n (sparse)
    norm2x_squared = np.array(norm2x_squared).ravel()
    norm2x = np.sqrt(norm2x_squared)
    norm2x_safe = norm2x + 1e-16

    # Normalization of the X columns
    inv_norms = 1.0 / norm2x_safe
    D = diags(inv_norms)
    Xn = X @ D
    # Solve ||X-ZSHO||_F with HOHO^T = D
    HO = orthNNLS(X, ZS, Xn)
    # Transposition
    Z = HO.T
    w, v = extract_w_v(Z)
    return w, v, Z


def initialize_w_values(Xl, v):
    n = len(v)
    w = np.zeros(n)
    # compute the degrees of the nodes
    degrees = np.array(Xl.sum(axis=1)).ravel()
    # compute the sum of node degrees for each community
    d_r = np.bincount(v, weights=degrees)
    w = degrees / d_r[v]
    return w

def extract_w_v(W):
    """ Extracts w and v from W."""
    w = np.max(W, axis=1)
    r = W.shape[1]
    v = np.argmax(W, axis=1)

    # attribuer aléatoirement un entier entre 0 et r-1 pour ces lignes
    rows_all_zero = np.all(W == 0, axis=1)
    v[rows_all_zero] = np.random.randint(0, r, size=np.sum(rows_all_zero))
    return w, v
def orthNNLS(M, U, Mn=None):
    """
    Solves the following optimization problem:
    min_{norm2v >= 0, V >= 0 and VV^T = D} ||M - U * V||_F^2

    Parameters:
        M (numpy.ndarray or csr_matrix): Matrix M of size (m, n).
        U (numpy.ndarray ): Matrix U of size (m, r).
        Mn (numpy.ndarray or csr_matrix, optional): Normalized columns of M. If None, it will be computed.

    Returns:
        V (numpy.ndarray): The matrix V of size (r, n) that approximates M.
        norm2v (numpy.ndarray): The squared norms of the columns of V.
    see F. Pompili, N. Gillis, P.-A. Absil and F. Glineur, "Two Algorithms for
    Orthogonal Nonnegative Matrix Factorization with Application to
    Clustering", Neurocomputing 141, pp. 15-25, 2014.
    """

    if Mn is None:
        # Normalize columns of M
        if issparse(M):
            norm2x_squared = M.multiply(M).sum(axis=0)  # matrice 1 x n (sparse)
            norm2x_squared = np.array(norm2x_squared).ravel()
            norm2x = np.sqrt(norm2x_squared)
            norm2x_safe = norm2x + 1e-16

            # Créer matrice diagonale inverses des normes
            inv_norms = 1.0 / norm2x_safe
            D = diags(inv_norms)  # matrice diagonale sparse (n x n)

            # Normaliser X par colonnes : multiplication à droite
            Mn = M @ D
        else:
            norm2m = np.sqrt(np.sum(M ** 2, axis=0))  # norm2m is the L2 norm of each column of M
            Mn = M * (1 / (norm2m + 1e-16))  # Avoid division by zero

    m, n = Mn.shape
    m_, r = U.shape

    # Normalize columns of U
    norm2u = np.sqrt(np.sum(U ** 2, axis=0))  # norm2u is the L2 norm of each column of U
    Un = U * (1 / (norm2u + 1e-16))  # Avoid division by zero
    if issparse(M):
        Mn = Mn.tocsc()
        M = M.tocsc()

    # Calculate the matrix A, which is the angles between columns of M and U
    A = Mn.T @ Un  # A is (n, r), matrix of angles

    # Find the index of the maximum value in each row of A (best column of U to approximate each column of M)
    b = np.argmax(A, axis=1)  # Indices of the best matching column in U

    # Initialize V with zeros
    V = np.zeros((r, n))

    # Assign the optimal weights to V(b(i), i)
    for i in range(n):
        if issparse(M):
            V[b[i], i] = (M[:, i].T @ U[:, b[i]] )[0] / norm2u[b[i]] ** 2
        else:
            V[b[i], i] = np.dot(M[:, i].T, U[:, b[i]]) / norm2u[b[i]] ** 2

    return V


# ------------------------------------------------
# -------------- UTILS ------------------
# ------------------------------------------------
def compute_error(normX, S):
    """ Computes error ||X - WSW'||_F """
    error = np.sqrt(1e-9 + normX ** 2 - np.linalg.norm(S, 'fro') ** 2) 
    return error
def compute_error_stupid(X,w,v,S):
    """ Computes error ||X - WSW'||_F """
    r = S.shape[0]
    n = X.shape[0]
    W = np.zeros((n, r))
    for i in range(n):
        W[i, v[i]] = w[i]

    error = np.linalg.norm(X-W@S@(W.T), 'fro')
    return error

def PowMethOTRISYMNMFFixed(X, r, v, w, maxiter, timelimit):
    t0 = time.process_time()
    e, t = [], []
    iter = 1
    n = X.shape[0]

    # Normalization of W
    nw = np.zeros(r)
    for i in range(n):
        nw[v[i]] += w[i] ** 2
    nw = np.sqrt(nw)

    denom = nw[v]
    mask = denom != 0
    w[mask] /= denom[mask]
    w[~mask] = 0

    # Main loop
    while iter <= maxiter and (time.process_time() - t0) <= timelimit:
        w_prev = [vec.copy() for vec in w]

        for i in range(1, r + 1):
            Ii = np.where(v == i)[0]
            x = np.zeros(len(Ii))

            for j in range(1, r + 1):
                Ij = np.where(v == j)[0]
                if len(Ii) > 0 and len(Ij) > 0:
                    vec = (X[np.ix_(Ii, Ij)] @ w[Ij]).ravel()
                    term = (w[Ii].T @ X[np.ix_(Ii, Ij)] @ w[Ij])
                    x = x + term * vec  # ici x et vec ont la même taille

            if np.linalg.norm(x) != 0:
                w[Ii] = x / np.linalg.norm(x)
        if sum(sum(abs(w_prev - w) / np.linalg.norm(w, 'fro'))) < 1e-7:
            break

        t.append(time.process_time() - t0)

        iter += 1

    # Normalization of W

    nw = np.zeros(r)
    for i in range(n):
        nw[v[i]] += w[i] ** 2
    nw = np.sqrt(nw)

    denom = nw[v]
    mask = denom != 0
    w[mask] /= denom[mask]
    w[~mask] = 0

    return w, v
