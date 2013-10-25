import sys
import os.path as op
from toolshed import reader
from matplotlib import pyplot as plt
from mpltools import style
style.use('ggplot')
from glob import glob
import operator
import numpy as np


path = "%s/ns_{ns}/rho_{rho}/w_{w}/" % (sys.argv[1])

COLORS = '#348ABD #7A68A6 #A60628 #467821 #CF4457 #188487 #E24A33'.split()



f, axes = plt.subplots(nrows=4, ncols=2, figsize=(10, 5))

def count_lt(fname, n_probes=2, check=operator.eq, p_cutoff=1e-5):
    return sum(1 for d in reader(fname) if float(d['p']) <= p_cutoff
                                           and check(int(d['n_probes']),
                                                   n_probes))

def basename(f):
    return op.basename(f).rstrip('.sim.covs.pvals.bed').rsplit('-', 1)[0]


regions = (2, 3, 4, 5)

rho = 0.3
for ix, region_size in enumerate(regions):
    check = operator.ge if region_size == 5 else operator.eq
    check = operator.eq
    for ins, ns in enumerate((10, 20, 40)):

        for iw, w in enumerate((0, 0.8)):
            directory = path.format(**locals())
            beds = sorted([x for x in glob("%s/*.bed" % directory)
                                                    if not ("sds_3" in x or
                                                            "skat" in x or
                                                            "bump" in x)])
            counts = [(f, count_lt(f, region_size, check)) for f in beds]

            ax = axes[ix, iw]
            leftish = len(counts) + 1.0

            print (ix, iw), ins, ins + np.arange(len(counts)) / leftish, [c[1]
                    for c in counts]

            rects = ax.bar(ins + np.arange(len(counts), dtype='f') / leftish,
                    [c[1] for c in counts],
                    width=1.0 / leftish - 0.015)
            for ir, r in enumerate(rects):
                r.set_facecolor(COLORS[ir])


rows, cols = axes.shape
for i in range(rows):
    for j in range(cols):
        axes[i, j].set_xticks([])
        if i == rows - 1:
            axes[i, j].set_xticks((0.5, 1.5, 2.5))
            axes[i, j].set_xticklabels(('10 samples', '20 samples', '40 samples'))
        if j == 0:
            axes[i, j].set_ylabel('%i probes' % regions[i])
            axes[i, j].set_ylim((0, 50))
        if j == 1:
            axes[i, j].set_ylim((0, 9000))

        if i == 0:
            axes[i, j].set_title('true +' if j == 1 else 'false +')


axes[0, 0].legend(rects, [basename(b) for b in beds], mode="expand", ncol=2)
plt.show()
