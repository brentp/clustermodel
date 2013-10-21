# from: http://nbviewer.ipython.org/urls/raw.github.com/EnricoGiampieri/dataplot/master/statplot.ipynb
from scipy.stats import gaussian_kde, sem
import numpy as np
import pylab as plt
import sys

def half_horizontal_bar(data, pos, left=False, dmin=0, dmax=1, **kwargs):
    n_bins = 40
    ax = kwargs.pop('ax', plt.gca())

    bins = np.linspace(dmin, dmax, n_bins + 1)
    #bins = np.linspace(data.min() - 0.01, data.max() + 0.01, n_bins + 1)

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
    dmax = min(1, dmax)

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
