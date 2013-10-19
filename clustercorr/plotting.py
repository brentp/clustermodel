# from: http://nbviewer.ipython.org/urls/raw.github.com/EnricoGiampieri/dataplot/master/statplot.ipynb
from scipy.stats import gaussian_kde,sem
import numpy as np
import pylab as plt

def half_horizontal_bar(data, pos, left=False, **kwargs):
    n_bins = 40
    ax = kwargs.pop('ax', plt.gca())

    bins = np.linspace(data.min() - 0.01, data.max() + 0.01, n_bins + 1),
    bins = np.linspace(0, 1, n_bins + 1)

    counts, edges = np.histogram(data, bins=bins, density=True)

    bsize = edges[1] - edges[0]
    counts /= float(counts.sum())
    if left:
        counts *= -1

    pos += (-0.0002 if left else 0.0002)
    ax.barh(edges[:n_bins], counts, bsize, left=pos, **kwargs)

    #return ax.fill_betweenx(x,pos,pos+v,**kwargs)

def hbar_plot(data1,classes=None,data2=None,**kwargs):
    ax = kwargs.get('ax', plt.gca())
    positions=range(len(data1))
    data2 = data2 if data2 is not None else data1
    classes = classes if classes is not None else positions
    assert len(classes)==len(data1) and len(classes)==len(data2)
    for pos,key in zip(positions, classes):
        try:
            d1,d2=data1[key], data2[key]
        except TypeError:
            d1,d2=data1[pos], data2[pos]
        color1=kwargs.pop('color1','b')
        color2=kwargs.pop('color2','b' if data1 is data2 else 'r')
        half_horizontal_bar(d1, pos, False, facecolor=color1)
        half_horizontal_bar(d2, pos, True, facecolor=color2)

    ax.set_xticks(positions)
    ax.set_xticklabels([str(i) for i in classes])
