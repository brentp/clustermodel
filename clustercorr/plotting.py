# from: http://nbviewer.ipython.org/urls/raw.github.com/EnricoGiampieri/dataplot/master/statplot.ipynb
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
from collections import OrderedDict

def half_horizontal_bar(data, pos, left=False, dmin=0, dmax=1, **kwargs):
    n_bins = 40 if len(data) > 15 else 20
    ax = kwargs.pop('ax', plt.gca())
    bins = np.linspace(dmin, dmax, n_bins + 1)

    counts, edges = np.histogram(data, bins=bins, density=True)
    counts = (0 + counts)

    bsize = edges[1] - edges[0]
    counts /= (2.5 * float(counts.max()))

    if left:
        counts *= -1

    pos += (-0.0002 if left else 0.0002)
    return ax.barh(edges[:n_bins], counts, bsize, left=pos, **kwargs)[0]

COLORS = '#348ABD #7A68A6 #A60628 #467821 #CF4457 #188487 #E24A33'.split()

def hbar_plot(data1, classes=None, data2=None, chrom='', **kwargs):
    ax = kwargs.get('ax', plt.gca())
    positions = range(len(data1))
    data2 = data2 if data2 is not None else data1
    classes = classes if classes is not None else positions
    assert len(classes) == len(data1) and len(classes) == len(data2)

    dmin = min(data1[key].min() for key in classes)
    dmax = max(data1[key].max() for key in classes)

    dmin = min(dmin, min(data2[key].min() for key in classes)) - 0.1
    dmax = max(dmax, min(data2[key].max() for key in classes)) + 0.1

    #dmin = max(0, dmin)
    #dmax = min(1, dmax)

    for pos, key in zip(positions, classes):
        try:
            d1, d2 = data1[key], data2[key]
        except TypeError:
            d1, d2 = data1[pos], data2[pos]
        shape1 = half_horizontal_bar(d1, pos, True, facecolor=COLORS[0], dmin=dmin,
                dmax=dmax)
        shape2 = half_horizontal_bar(d2, pos, False, facecolor=COLORS[3], dmin=dmin,
                dmax=dmax)

    ax.set_ylim(dmin, dmax)
    ax.set_xticks(positions)
    ax.set_xlim(-0.5, max(positions) + 0.5)
    if chrom: chrom += ":"
    if isinstance(classes[0], int):
        lbls = ["%s%s" % (chrom, "{:,}".format(i)) for i in classes]
    else:
        lbls = [str(s) for s in classes]
    ax.set_xticklabels(lbls,
            rotation=15 if len(classes) > 8 else 0)
    return shape1, shape2

def plot_hbar(covs, cluster_df, covariate, chrom, res, png):
    from matplotlib import pyplot as plt
    group = getattr(covs, covariate)
    grps = sorted(list(set(group)))

    ax = plt.gca()
    cdf = cluster_df.T
    #cdf = 1 / (1 + np.exp(-cdf))

    d1 = OrderedDict(cdf.ix[group == grps[0], :].T.iterrows())
    d2 = OrderedDict(cdf.ix[group == grps[1], :].T.iterrows())

    r1, r2 = hbar_plot(d1, list(cdf.columns), d2, ax=ax, chrom=chrom)

    ax.legend((r1, r2),
             ("%s - %s" % (covariate, grps[0]),
              "%s - %s" % (covariate, grps[1])),
              loc='upper left')

def plot_continuous(covs, cluster_df, covariate, chrom, res, png):
    from matplotlib import pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(ncols=cluster_df.shape[0])
    cdf = cluster_df.T
    #cdf.columns = ['%s:%s' % (chrom, "{:,}".format(p)) for p in cdf.columns]

    for i, c in enumerate(cdf.columns):
        ax = axes[i]
        cname = "%s:%s" % (chrom, "{:,}".format(c))
        ax.plot(cdf[c], getattr(covs, covariate), marker='o', ls='none')
        ax.set_title(cname)
        ax.set_ylabel(covariate)
        ax.set_xlabel('methylation')
    return fig


def plot_dmr(covs, cluster_df, covariate, chrom, res, png):
    from matplotlib import pyplot as plt
    import numpy as np
    from pandas.tools.plotting import parallel_coordinates

    cdf = cluster_df.T
    cdf.columns = ['%s:%s' % (chrom, "{:,}".format(p)) for p in cdf.columns]
    cdf = 1 / (1 + np.exp(-cdf))
    cdf['group'] = getattr(covs, covariate)

    ax = plt.gca()

    if cdf.group.dtype == float:
        ax = parallel_coordinates(cdf, 'group', ax=ax)
        ax.get_legend().set_visible(False)
    else:
        ax = parallel_coordinates(cdf, 'group', colors=('#764AE7', '#E81C0E'),
                ax=ax)
        lbls = ax.get_legend().get_texts()

        for lbl in lbls:
            lbl.set_text(covariate + ' ' + lbl.get_text())

    if len(cdf.columns) > 6:
        ax.set_xticklabels([x.get_text() for x in ax.get_xticklabels()],
                          rotation=10)

    ax.set_ylabel('methylation')
