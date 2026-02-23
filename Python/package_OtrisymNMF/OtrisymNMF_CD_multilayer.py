import numpy as np
import time
from scipy.sparse import issparse, find
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import normalized_mutual_info_score
from .SVCA import SVCA
from .orthNNLS import orthNNLS
import math
from package_OtrisymNMF.aggregation_way.co_association import clustering_coassociation, MatrixFreeSpectralClusteringCoAssociation, \
    USENC_ConsensusFunction
from package_OtrisymNMF.aggregation_way.QR_embeddings import clustering_QR_embeddings, nmf_consensus,consensus_kmeans

import numpy as np
# import Cluster_Ensembles as CE


# -----------------------------------------------
# ---------------- OtrisymNMF_CD ----------------
# -----------------------------------------------
def OtrisymNMF_CD(X, r, numTrials=3, maxiter=1000, delta=1e-7, time_limit=300, init_method='Scalable_co', verbosity=0,
                  init_seed=None, true_labels = None):
    errors = []

    start_time = time.time()


    # dense format numpy
    if issparse(X):
        X = X.toarray()

    L = X.shape[0]
    n = X.shape[1]
    error_best = float('inf')

    w_best = np.zeros((L, n))
    S_best = np.zeros((L, r, r))
    v_best = np.zeros(n)

    # Precomputations
    normX = np.zeros(L)
    for l in range(L):
        normX[l] = np.linalg.norm(X[l], 'fro')

    if verbosity > 0:
        print(f'Running {numTrials} Trials in Series')

    for trial in range(numTrials):
        if init_method is None:
            init_algo = 'onelayer'
        else:
            init_algo = init_method

        if init_seed is not None:
            init_seed += 10 * trial
            np.random.seed(init_seed)

        if init_algo == 'random':
            base = np.arange(r)
            rest = np.random.randint(0, r, size=n - r)
            v = np.concatenate([base, rest])
            np.random.shuffle(v)
            w = np.zeros((L, n))
            for l in range(L):
                w[l] = initialize_w_values(X[l], v)
            w, v = PowMethOTRISYMNMFFixed(X, r, v, w, maxiter=200, timelimit=200)

        elif init_algo == 'onelayer':
            w = np.zeros((L, n))
            lr = np.random.randint(0, L) # Choose a random layer
            w[lr], v, _ = initialize_W_onelayer(X[lr], r)
            # default value degree of the node / sum degrees same community
            for l in range(L):
                if l == lr:
                    continue
                w[l] = initialize_w_values(X[l], v)
            w, v = PowMethOTRISYMNMFFixed(X, r, v, w, maxiter=200, timelimit=200)


        else:

            w, v = initialize_W_alllayers(X, r, init_method)


        # Normalization of W
        for l in range(L):
            nw = np.zeros(r)
            for i in range(n):
                nw[v[i]] += w[l][i] ** 2
            nw = np.sqrt(nw)

            denom = nw[v]
            mask = denom != 0
            w[l][mask] /= denom[mask]
            w[l][~mask] = 0
        # Compute of S
        S = update_S(X, r, w, v)

        if verbosity and (true_labels is not None):
            print('NMI initialisation : ', normalized_mutual_info_score(true_labels, v))
            print('Time',time.time()-start_time)
        prev_error = 0
        for l in range(L):
            prev_error += compute_error(normX[l], S[l])
        error = prev_error
        init_error = (error / L)




        errors.append(init_error)
        for iteration in range(maxiter):
            if time.time() - start_time > time_limit:
                print('Time limit passed')
                break

            w, v = update_W(X, S, w, v)

            S = update_S(X, r, w, v)

            prev_error = error
            error = 0
            for l in range(L):

                error += compute_error(normX[l], S[l])


            if error < delta or abs(prev_error - error) < delta:
                break

        if error < error_best:
            w_best, v_best, S_best, error_best = (
                w.copy(), v.copy(), S.copy(), error
            )

        if verbosity > 0:
            print(f'Trial {trial + 1}/{numTrials} with {init_algo}: Error {error:.4e} | Best: {error_best:.4e}')
            if true_labels is not None:
                print('NMI : ', normalized_mutual_info_score(true_labels, v_best))
            print('Time', time.time()-start_time)

    return w_best, v_best, S_best, errors

# ------------------------------------------------
# -------------- Initialization ------------------
# ------------------------------------------------
def initialize_W_onelayer(Xl,r):

    options = {'average': 1}
    n = Xl.shape[0]
    p = max(2, math.floor(0.1 * n / r))
    WO, K = SVCA(Xl, r, p, options=options)
    norm2x = np.sqrt(np.sum(Xl ** 2, axis=0))
    S_sum_n = Xl * (1 / (norm2x + 1e-16))
    HO = orthNNLS(Xl, WO, S_sum_n)
    W = HO.T
    w, v = extract_w_v(W)
    return w, v, W

def initialize_w_values(Xl,v):
    n = len(v)
    w = np.zeros(n)
    # compute the degrees of the nodes
    degrees = Xl.sum(axis=1)
    # compute the sum of node degrees for each community
    d_r = np.bincount(v, weights=degrees)
    w = degrees / d_r[v]
    return w

def initialize_W_alllayers(X, r, init_method):
    n = X.shape[1]
    L = X.shape[0]


    W = np.zeros((L, n, r))
    w = np.zeros((L, n))
    v_init = np.zeros((L, n))
    # Find w and v for each layer
    for l in range(L):
        w[l], v_init[l], W[l] = initialize_W_onelayer(X[l], r)

    if init_method == 'co-association':
        C = clustering_coassociation(W) # calcul de la matrice de co-association à partir de W
        model = SpectralClustering(n_clusters=r, affinity='precomputed', random_state=42)
        labels_final = model.fit_predict(C)

    elif init_method == 'QR':
        C = clustering_QR_embeddings(W)
        clustering_model = KMeans(n_clusters=r, n_init=50, random_state=42)
        labels_final = clustering_model.fit_predict(C)

    elif init_method =='NMF':
        labels_final= nmf_consensus(v_init, r)
    elif init_method =='simple':
        labels_final=consensus_kmeans(v_init, r, normalize_rows=True)
    elif init_method =='Scalable_co':
        method = MatrixFreeSpectralClusteringCoAssociation(v_init, r)
        labels_final = method.fit_predict()

    elif init_method == 'USENC':
        labels_final = USENC_ConsensusFunction(v_init.T, r)

    for l in range(L):
        w[l] = initialize_w_values(X[l], labels_final)

    #w_new, v = PowMethOTRISYMNMFFixed(X, r, labels_final, w, maxiter=200, timelimit=200)

    return w, labels_final


def extract_w_v(W):
    """ Extracts w and v from W."""
    w = np.max(W, axis=1)
    r = W.shape[1]
    v = np.argmax(W, axis=1)

    # attribuer aléatoirement un entier entre 0 et r-1 pour ces lignes
    rows_all_zero = np.all(W == 0, axis=1)
    v[rows_all_zero] = np.random.randint(0, r, size=np.sum(rows_all_zero))
    return w, v


def update_W(X, S, w, v):
    n = X.shape[1]
    r = S.shape[1]
    L = X.shape[0]
    """
    Pré-calculs pour éviter une double boucle sur "n"
    wp2 n'est plus un vecteur, mais une matrice puisque w est propre à chaque couche
    """
    wp2 = np.zeros((L, r))

    for l in range(L):
        for k in range(r):
            for p in range(n):
                wp2[l, k] += (w[l, p] * S[l, v[p], k]) ** 2

    # Mise à jour de W
    for i in range(n):
        vi_new = -1
        wi_new = -1 * np.ones(L)
        wi = np.zeros((L, n))
        f_new = np.inf

        Xii = np.zeros(L)
        for l in range(L):
            Xii[l] = X[l][i][i]

        for k in range(r):
            wi = w.copy()
            vi = v.copy()
            vi[i] = k
            erreur = 0
            # Boucle sur chaque couche
            for l in range(L):
                c3 = S[l, k, k] ** 2
                c1 = 2 * (wp2[l][k] - (w[l][i] * S[l, v[i], k]) ** 2) - 2 * S[l][k][k] * Xii[l]

                c0 = -4 * sum(X[l][i][p] * w[l][p] * S[l][v[p]][k] for p in np.nonzero(X[l][i, :])[0] if p != i)

                # Résolution des racines avec la méthode de Cardan
                roots = cardan_depressed(4 * c3, 2 * c1, c0)

                # Trouver la meilleure solution positive pour w_l(i, k_i)
                x = np.sqrt(r / n)
                min_value = c3 * (x ** 4) + c1 * (x ** 2) + c0 * x
                for sol in roots:
                    value = c3 * (sol ** 4) + c1 * (sol ** 2) + c0 * sol
                    if sol > 0 and value < min_value:
                        x, min_value = sol, value

                wi[l, i] = x

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
            w[l][i] = wi_new[l]


    # Normalization of W
    for l in range(L):
        nw = np.zeros(r)
        for i in range(n):
            nw[v[i]] += w[l][i] ** 2
        nw = np.sqrt(nw)

        denom = nw[v]
        mask = denom != 0
        w[l][mask] /= denom[mask]
        w[l][~mask] = 0

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


def update_S(X, r, w, v):
    L = X.shape[0]
    S = np.zeros((L, r, r))
    for l in range(L):

        # Get the row indices, column indices, and values of the non-zero elements in the sparse matrix X
        i, j, val = find(X[l])

        # Loop through the non-zero elements of X
        for k in range(len(val)):
            S[l, v[i[k]], v[j[k]]] += w[l, i[k]] * w[l, j[k]] * val[k]
    return S


def compute_error(normX, S):
    """ Computes error ||X - WSW'||_F / ||X||_F."""
    error = np.sqrt(1e-9 + normX ** 2 - np.linalg.norm(S, 'fro') ** 2) / normX
    return error



def PowMethOTRISYMNMFFixed(X, r, v, w, maxiter, timelimit):
    t0 = time.process_time()
    e, t = [], []
    iter = 1
    L = X.shape[0]
    n = X.shape[1]

    # Normalization of W
    for l in range(L):
        nw = np.zeros(r)
        for i in range(n):
            nw[v[i]] += w[l][i] ** 2
        nw = np.sqrt(nw)

        denom = nw[v]
        mask = denom != 0
        w[l][mask] /= denom[mask]
        w[l][~mask] = 0

    # Main loop
    while iter <= maxiter and (time.process_time() - t0) <= timelimit:
        w_prev = w.copy()
        for l in range(X.shape[0]):

            for i in range(1, r+1):
                Ii = np.where(v == i)[0]
                x = np.zeros(len(Ii))
                Xl = X[l]

                for j in range(1, r+1):
                    Ij = np.where(v == j)[0]
                    if len(Ii) > 0 and len(Ij) > 0:
                        vec = (Xl[np.ix_(Ii, Ij)] @ w[l, Ij]).ravel()
                        term = (w[l, Ii].T @ Xl[np.ix_(Ii, Ij)] @ w[l, Ij])
                        x = x + term * vec  # ici x et vec ont la même taille

                if np.linalg.norm(x) != 0:
                    w[l, Ii] = x / np.linalg.norm(x)
        if sum(sum(abs(w_prev-w)/np.linalg.norm(w, 'fro'))) < 1e-7:
            break

        t.append(time.process_time() - t0)


        iter += 1

    # Normalization of W
    for l in range(L):
        nw = np.zeros(r)
        for i in range(n):
            nw[v[i]] += w[l][i] ** 2
        nw = np.sqrt(nw)

        denom = nw[v]
        mask = denom != 0
        w[l][mask] /= denom[mask]
        w[l][~mask] = 0

    return w, v
