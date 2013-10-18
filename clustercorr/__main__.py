import sys
import pandas as pd
from aclust import aclust
from clustercorr import feature_gen, cluster_to_dataframe, clustered_model

def clustermodel(fcovs, fmeth, model, max_dist=500, linkage='complete',
        rho_min=0.3, min_clust_size=2, sep="\t",
        liptak=False, bumping=False, gee_args=(), skat=False):
    # an iterable of feature objects
    feature_iter = feature_gen(fmeth, rho_min=rho_min)

    assert min_clust_size > 1

    cluster_gen = (c for c in aclust(feature_iter, max_dist=max_dist,
                                     max_skip=1, linkage=linkage)
                    if len(c) >= min_clust_size)

    covs = pd.read_table(fcovs, index_col=0, sep=sep)

    for cluster in cluster_gen:

        # we turn the cluster list into a pandas dataframe with columns
        # of samples and rows of probes. these must match our covariates
        cluster_df = cluster_to_dataframe(cluster, columns=covs.index)

        # now we want to test a model on our clustered dataset.
        res = clustered_model(covs, cluster_df, model, gee_args=gee_args,
                liptak=liptak, bumping=bumping, skat=skat)
        res['chrom'] = cluster[0].group
        res['start'] = cluster[0].start
        res['end'] = cluster[-1].end
        res['n_probes'] = len(cluster)
        yield res

def main_example():
    fcovs = "clustercorr/tests/example-covariates.txt"
    fmeth = "clustercorr/tests/example-methylation.txt.gz"
    model = "methylation ~ disease + gender"

    for cluster_p in clustermodel(fcovs, fmeth, model):
        if cluster_p['p'] < 1e-4:
            print cluster_p

def main(args=sys.argv[1:]):
    import argparse
    p = argparse.ArgumentParser(__doc__)
    mp = p.add_argument_group('modeling choices (choose one or specify a '
            'mixed-model using lme4 syntax)')
    group = mp.add_mutually_exclusive_group()

    group.add_argument('--skat', action='store_true')
    group.add_argument('--gee-args',
                       help='comma-delimited correlation-structure, variable')
    group.add_argument('--liptak', action="store_true")
    group.add_argument('--bumping', action="store_true")


    cp = p.add_argument_group('clustering parameters')
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

    p.add_argument('model',
                   help="model in R syntax, e.g. 'methylation ~ disease'")
    p.add_argument('covs', help="tab-delimited file of covariates: shape is "
                   "n_samples * n_covariates")
    p.add_argument('methylation', help="tab-delimited file of methylation"
                   " rows of this file must match the columns of `covs`"
                   " shape is n_probes * n_samples")
    a = p.parse_args(args)

    if a.gee_args is not None:
        method = 'gee:' + a.gee_args
        a.gee_args = a.gee_args.split(",")
    else:
        if a.liptak: method = 'liptak'
        elif a.bumping: method = 'bumping'
        elif a.skat: method = 'skat'
        else:
            assert "|" in a.model
            method = "mixed-model"

    fmt = "{chrom}\t{start}\t{end}\t{coef}\t{p}\t{n_probes}\t{model}\t{method}"
    print "#" + fmt.replace("}", "").replace("{", "")
    for c in clustermodel(a.covs, a.methylation, a.model,
                          max_dist=a.max_dist,
                          linkage=a.linkage,
                          rho_min=a.rho_min,
                          min_clust_size=a.min_cluster_size,
                          liptak=a.liptak,
                          bumping=a.bumping,
                          gee_args=a.gee_args,
                          skat=a.skat):
        c['method'] = method
        print fmt.format(**c)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "example":
        sys.exit(main_example())

    main()
