import numpy as np
import pandas as pd
import scipy.stats as ss
from aclust import aclust

from . import feature

choice = np.random.choice

def generate(clust_gen, samples, n=20, w=1):
    """
    if w == 0, then this returns random data
    """
    samples = np.asarray(samples)
    for cluster in clust_gen:
        N = len(cluster[0].values)
        assert 2 * n < N

        n_probes = len(cluster)
        new_data = np.zeros((n_probes, 2 * n))

        # choose a ranomd probe from the set.
        # the values across the cluster will be determined
        # by the values in this randomly chose probe.
        i = choice(range(n_probes))[0]
        c = cluster[i]

        idxs = np.arange(N)

        # just pull based on the index. so we need to sort the values
        # as well.
        idx_order = np.argsort(c.values)

        ords = np.arange(1, N + 1) / (N + 1.0)
        ords = (1.0 - ords)**w
        h_idxs = choice(idxs, replace=False, p=ords/ords.sum(), size=n)

        idxs = np.setdiff1d(idxs, h_idxs, assume_unique=True)
        idxs.sort()

        ords = np.arange(1, N + 1 - n) / (N + 1.0 - n)
        assert ords.shape == idxs.shape
        ords = (ords)**w
        l_idxs = choice(idxs, replace=False, p=ords/ords.sum(), size=n)

        assert len(np.intersect1d(h_idxs, l_idxs)) == 0

        for j in range(n_probes):
            # sort then pull
            new_data[j,:n] = cluster[j].values[idx_order][h_idxs]
            new_data[j,n:] = cluster[j].values[idx_order][l_idxs]

        new_data = pd.DataFrame(new_data)
        new_data.index = ["%s:%i" % (c.group, c.end) for c in cluster]
        new_data.columns = list(samples[idx_order][h_idxs]) + \
                           list(samples[idx_order][l_idxs])
        yield new_data

def gen_clusters_from_data(fmeth, min_clust_size=3, max_clust_size=200, n=20, w=0.5):
    header = open(fmeth).readline().rstrip().split("\t")[1:]
    clust_gen = (c for c in aclust(feature.feature_gen(fmeth),
                                   max_dist=500,
                                   max_skip=1,
                                   linkage=0.5)
                              if len(c) >= min_clust_size
                             and len(c) <= max_clust_size)
    for cluster in generate(clust_gen, header, n=n, w=w):
        yield cluster


def gen_arma(n_probes=20000, n_patients=80, corr=0.2, df=2, scale=0.2):
    # these parameters are taken from the Bump Hunting Paper
    from statsmodels.tsa.arima_process import ArmaProcess
    sigma = 0.5
    rvs = ss.norm(df, loc=0.05 / sigma, scale=scale).rvs
    corr = -abs(corr)
    return np.column_stack([
        sigma * ArmaProcess([1, corr], [1])
                           .generate_sample(n_probes=n_probes, distrvs=rvs)
                           for i in range(n_patients)])

if __name__ == "__main__":

    fmeth = 'work/chr1.M.txt'
    n = 20

    for cluster in gen_clusters_from_data(fmeth, n=n, min_clust_size=3, w=0.9):
        print cluster.ix[:, :n].mean(axis=1)
        print cluster.ix[:, n:].mean(axis=1)
        print cluster.ix[:, :n].mean(axis=1) - cluster.ix[:, n:].mean(axis=1)
        print len(cluster)
        print

