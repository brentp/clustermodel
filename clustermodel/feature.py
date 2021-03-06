import sys
import numpy as np
import scipy.stats as ss
import pandas as pd
from toolshed import reader

class ClusterFeature(object):
    __slots__ = "group start end values rho_min weights".split()

    def __init__(self, group, start, end, values, rho_min=0.25, weights=None):
        self.group, self.start, self.end, self.values = group, start, end, values
        self.weights = weights
        self.rho_min = rho_min

    def distance(self, other):
        if self.group != other.group: return sys.maxint
        if other.start > self.end:
            return other.start - self.end
        elif self.start > other.end:
            return self.start - other.end
        return 0

    def is_correlated(self, other):
        rho, p = ss.spearmanr(self.values, other.values)
        return rho > self.rho_min

    def __repr__(self):
        c = self.__class__.__name__
        return "%s(%s:%s-%s [%i values])" % (c, self.group, self.start,
                                             self.end, len(self.values))

def row_handler(tokens):
    try:
        chrom, pos = tokens[0].split(":" if ":" in tokens[0] else "_")
    except ValueError:
        chrom, pos = "chrom", int(tokens[0])
    return (chrom, int(pos) - 1, int(pos), np.array([float(x or 'nan')
                                                     for x in tokens[1:]]))

def cluster_to_dataframe(cluster, columns=None, weights=False):
    if weights:
        df = pd.DataFrame([c.weights for c in cluster],
            index=["%s:%i" % (c.group, c.end) for c in cluster])
    else:
        df = pd.DataFrame([c.values for c in cluster],
            index=["%s:%i" % (c.group, c.end) for c in cluster])
    if columns is not None:
        df.columns = columns
    return df

def feature_gen(fname, row_handler=row_handler, feature_class=ClusterFeature, sep="\t",
        rho_min=0.3, skip_first_row=True, weights=None):
    """

    Parameters
    ----------
    fname : str
        file name containing methylation data

    row_handler: function
        function that takes a list of values for each line in `fname` and
        returns a tuple of chrom, start, end, values. e.g.
        def row_handler(tokens):
            chrom, pos = tokens[0].split(":")
            return (chrom, int(pos) - 1, int(pos), map(float, values[1:]))

    feature_class: class
        a class derived from `ClusterFeature` that accepts
        chrom, start, end, values and has those atributes and fulfills
        the requirements of aclust.aclust.

    rho_min: float
        the minimum spearman's r between 2 sets of values for them to be
        considered as correlated
    """
    if weights is not None:
        weights = reader(weights, header=False, sep=sep)
    for i, toks in enumerate(reader(fname, header=False, sep=sep)):
        if i == 0 and skip_first_row:
            if weights is not None: next(weights)
            continue
        vals = row_handler(toks)
        if weights is not None:
            chrom, start, end, weight_vals = row_handler(next(weights))
            assert chrom == vals[0]
            assert start == vals[1], (vals[1], start)
        else:
            weight_vals = None
        yield feature_class(*vals, **{'rho_min': rho_min, 'weights':weight_vals} )
