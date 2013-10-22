# from: http://nbviewer.ipython.org/urls/raw.github.com/EnricoGiampieri/dataplot/master/statplot.ipynb
from scipy.stats import gaussian_kde, sem
import numpy as np
import pylab as plt
import sys

def half_horizontal_bar(data, pos, left=False, dmin=0, dmax=1, **kwargs):
    n_bins = 40
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
    return

    keep = ~np.isnan(data)
    x = bins
    v = gaussian_kde(data[keep]).evaluate(x)
    v = 0.33 * v/v.max() * (-1 if left else 1)

    kwargs['edgecolor'] = 'k'
    kwargs['alpha'] = 0.33
    kwargs['zorder'] = -1
    ax.fill_betweenx(x, pos, pos+v, **kwargs)

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

    dmin = max(0, dmin)
    #dmax = min(1, dmax)

    for pos, key in zip(positions, classes):
        try:
            d1, d2 = data1[key], data2[key]
        except TypeError:
            d1, d2 = data1[pos], data2[pos]
        color1 = kwargs.pop('color1', 'b')
        color2 = kwargs.pop('color2', 'b' if data1 is data2 else 'r')
        shape1 = half_horizontal_bar(d1, pos, False, facecolor=color1, dmin=dmin,
                dmax=dmax)
        shape2 = half_horizontal_bar(d2, pos, True, facecolor=color2, dmin=dmin,
                dmax=dmax)

    ax.set_ylim(dmin, dmax)
    ax.set_xticks(positions)
    ax.set_xlim(-0.5, max(positions) + 0.5)
    if chrom: chrom += ":"
    ax.set_xticklabels(["%s%s" % (chrom, "{:,}".format(i)) for i in classes],
            rotation=15 if len(classes) > 8 else 0)
    return shape1, shape2

def plot_hbar(covs, cluster_df, covariate, chrom, res, png):
    from matplotlib import pyplot as plt
    group = getattr(covs, covariate)
    grps = list(set(group))

    ax = plt.gca()
    cdf = cluster_df.T
    #cdf = 1 / (1 + np.exp(-cdf))

    d1 = dict(cdf.ix[group == grps[0], :].T.iterrows())
    d2 = dict(cdf.ix[group == grps[1], :].T.iterrows())

    r1, r2 = hbar_plot(d1, list(cdf.columns), d2, ax=ax, chrom=chrom)

    ax.legend((r1, r2),
             ("%s - %s" % (covariate, grps[0]),
              "%s - %s" % (covariate, grps[1])),
              loc='upper left')

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
