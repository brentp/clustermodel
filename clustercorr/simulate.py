import sys
from toolshed import nopen
import itertools

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
        yield simulate_cluster(cluster, samples, n, w)

def simulate_cluster(cluster, samples, n, w):
    N = len(cluster[0].values)
    assert 2 * n < N

    n_probes = len(cluster)
    new_data = np.zeros((n_probes, 2 * n))

    # choose a ranomd probe from the set.
    # the values across the cluster will be determined
    # by the values in this randomly chose probe.
    i = choice(range(n_probes))#[0] if n_probes > 1 else 0
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
    return new_data

def gen_clusters_from_data(fmeth, min_clust_size=3, max_clust_size=200, n=20, w=0.5):
    header = nopen(fmeth).readline().rstrip().split("\t")[1:]
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

def main(argv=sys.argv[1:]):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-samples", type=int, default=20,
            help="number of samples from each group to simulate")
    ap.add_argument("-w", type=float, default=0,
            help="weight parameter 0 is random, 1 is strong separation"
            " between simulated groups")
    ap.add_argument("methylation",
            help="input methylation data from which to simulate data")
    ap.add_argument("prefix",
            help="output prefix. prefix.meth.txt, prefix.covs.txt "
            "and prefix.clusters.txt" " be created")

    cp = ap.add_argument_group('clustering parameters')
    cp.add_argument('--rho-min', type=float, default=0.3,
                   help="minimum correlation to merge 2 probes")
    cp.add_argument('--min-cluster-size', type=int, default=2,
                    help="minimum cluster size on which to run model: "
                   "must be at least 2")
    cp.add_argument('--linkage', choices=['single', 'complete'],
                    default='complete', help="linkage method")

    cp.add_argument('--max-dist', default=500, type=int,
                    help="maximum distance beyond which a probe can not be"
                    " added to a cluster")

    args = ap.parse_args(argv)

    samples = np.array(nopen(args.methylation).readline().rstrip().split("\t")[1:])
    cluster_gen = (c for c in aclust(feature.feature_gen(args.methylation,
                                                         rho_min=args.rho_min),
                                     max_dist=args.max_dist,
                                     max_skip=1,
                                     linkage=args.linkage))
    def check_cluster(clust):
        return len(clust) >= args.min_cluster_size

    n, w = args.n_samples, args.w
    fh_meth = nopen("%s.meth.txt" % args.prefix, "w")
    fh_clst = nopen("%s.clusters.txt" % args.prefix, "w")
    with nopen("%s.covs.txt" % args.prefix, "w") as fh_covs:
        print >>fh_covs, "id\tcase"
        print >>fh_covs, "\n".join(
                         ["case_%i\t1" % i for i in range(n)] +
                         ["ctrl_%i\t0" % i for i in range(n)])
    print >>fh_meth, "probe\t" + "\t".join(
            ["case_%i" % i for i in range(n)] +
            ["ctrl_%i" % i for i in range(n)])

    print >>fh_clst, "chrom\tstart\tend\tw\tprobes"

    for i, ocluster in enumerate(cluster_gen):
        is_cluster = check_cluster(ocluster)
        cluster = simulate_cluster(ocluster, samples, n, w if is_cluster else 0)
        cluster.to_csv(fh_meth, sep="\t", index=True, header=False,
                float_format="%.3f")
        chrom, start, end = ocluster[0].group, ocluster[0].start, ocluster[-1].end
        probes = ",".join("%s:%i" % (o.group, o.end) for o in ocluster)
        if is_cluster:
            print >>fh_clst, "{chrom}\t{start}\t{end}\t{w}\t{probes}".format(**locals())

    print >>sys.stderr, "wrote:", fh_meth.name

if __name__ == "__main__":

    sys.exit(main())

    fmeth = 'work/chr1.M.txt'
    n = 20

    for cluster in gen_clusters_from_data(fmeth, n=n, min_clust_size=3, w=0.9):
        print cluster.ix[:, :n].mean(axis=1)
        print cluster.ix[:, n:].mean(axis=1)
        print cluster.ix[:, :n].mean(axis=1) - cluster.ix[:, n:].mean(axis=1)
        print len(cluster)
        print

