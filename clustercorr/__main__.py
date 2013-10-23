import sys
import numpy as np
import pandas as pd
import collections
from aclust import aclust
from .plotting import plot_dmr, plot_hbar, plot_continuous
from clustercorr import feature_gen, cluster_to_dataframe, clustered_model


def is_numeric(pd_series):
    if np.issubdtype(pd_series.dtype, int) or \
        np.issubdtype(pd_series.dtype, float):
        return len(pd_series.unique()) > 2
    return False

def print_compare(res, fmt, header=False):

    # lots of entries in res, we just want the ones with the
    # coefficients and p-values
    keys = [x for x in res.keys() if isinstance(res[x], tuple)
                                  and len(res[x]) == 2]
    cols = [x.replace(" ", "_") for x in keys]
    if header:
        h = "#" + fmt.replace("}", "").replace("{", "")
        h += "\t".join("p_{c}\tcoef_{c}".format(c=c) for c in cols)
        print h
    line = fmt.format(**res)
    vals = []
    for k in keys:
        vals.append("%.4g\t%.3f" % res[k])
    print line + "\t" + "\t".join(vals)

def clustermodel(fcovs, fmeth, model, max_dist=500, linkage='complete',
        rho_min=0.3, min_clust_size=2, sep="\t",
        outlier_sds=None,
        liptak=False, bumping=False, gee_args=(), skat=False, png_path=None,
        compare=False):
    # an iterable of feature objects
    feature_iter = feature_gen(fmeth, rho_min=rho_min)

    assert min_clust_size > 1

    cluster_gen = (c for c in aclust(feature_iter, max_dist=max_dist,
                                     max_skip=1, linkage=linkage)
                    if len(c) >= min_clust_size)

    covs = pd.read_table(fcovs, index_col=0, sep=sep)

    covariate = model.split("~")[1].split("+")[0].strip()

    if compare:
        compare = {"p": collections.defaultdict(list),
                   "coef": collections.defaultdict(list)}

    for cluster in cluster_gen:
        chrom = cluster[0].group

        # we turn the cluster list into a pandas dataframe with columns
        # of samples and rows of probes. these must match our covariates
        cluster_df = cluster_to_dataframe(cluster, columns=covs.index)

        # now we want to test a model on our clustered dataset.
        res = clustered_model(covs, cluster_df, model, gee_args=gee_args,
                liptak=liptak, bumping=bumping, skat=skat,
                outlier_sds=outlier_sds, compare=compare)
        res['chrom'] = cluster[0].group
        res['start'] = cluster[0].start
        res['end'] = cluster[-1].end
        res['n_probes'] = len(cluster)

        if compare:
            yield res
            continue


        if res['p'] < 1e-4 and png_path is not None:
            from matplotlib import pyplot as plt
            from mpltools import style
            style.use('ggplot')

            region = "{chrom}_{start}_{end}".format(**res)
            if png_path.endswith('show'):
                png = None
            elif png_path.endswith(('.png', '.pdf')):
                png = "%s.%s%s" % (png_path[:-4], region, png_path[-4:])
            elif png_path:
                png = "%s.%s.png" % (png_path.rstrip("."), region)

            if is_numeric(getattr(covs, covariate)):
                f = plot_continuous(covs, cluster_df, covariate, chrom, res, png)
            else:
                f = plt.figure(figsize=(11, 4))
                ax = f.add_subplot(1, 1, 1)
                if 'spaghetti' in png_path:
                    plot_dmr(covs, cluster_df, covariate, chrom, res, png)
                else:
                    plot_hbar(covs, cluster_df, covariate, chrom, res, png)
                plt.title('p-value: %.3g %s: %.3f' % (res['p'], covariate, res['coef']))
            f.set_tight_layout(True)
            if png:
                plt.savefig(png)
            else:
                plt.show()
        yield res


def main_example():
    fcovs = "clustercorr/tests/example-covariates.txt"
    fmeth = "clustercorr/tests/example-methylation.txt.gz"
    model = "methylation ~ disease + gender"

    for cluster_p in clustermodel(fcovs, fmeth, model):
        if cluster_p['p'] < 1e-5:
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
    group.add_argument('--compare', action="store_true",
                       help='run all methods and compare')

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
    p.add_argument('--png-path',
                   help="path to save a png of regions with low p-values. If "
                   "this ends with 'show', each plot will be shown in a window"
                   " if this contains the string 'spaghetti', it will draw a "
                   "a spaghetti plot, otherwise, it's a histogram plot")
    p.add_argument('--outlier-sds', type=float, default=30,
            help="remove points that are more than this many standard "
                 "deviations away from the mean (only usable with GEEs"
                 " and mixed-models which allow missing data")

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
        elif a.compare: method = 'compare'
        else:
            assert "|" in a.model
            method = "mixed-model"

    if a.compare:
        fmt = "{chrom}\t{start}\t{end}\t{n_probes}\t"
    else:
        fmt = "{chrom}\t{start}\t{end}\t{coef}\t{p}\t{n_probes}\t{model}\t{method}"
        print "#" + fmt.replace("}", "").replace("{", "")

    for i, c in enumerate(clustermodel(a.covs,
                          a.methylation, a.model,
                          max_dist=a.max_dist,
                          linkage=a.linkage,
                          rho_min=a.rho_min,
                          min_clust_size=a.min_cluster_size,
                          liptak=a.liptak,
                          bumping=a.bumping,
                          gee_args=a.gee_args,
                          skat=a.skat,
                          outlier_sds=a.outlier_sds,
                          png_path=a.png_path,
                          compare=a.compare)):
        if a.compare:
            print_compare(c, fmt, header=i==0)

        else:
            c['method'] = method
            print fmt.format(**c)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "example":
        sys.exit(main_example())
    if len(sys.argv) > 1 and sys.argv[1] == "simulate":
        from . import simulate
        sys.exit(simulate.main(sys.argv[2:]))

    main()
