from nose.tools import assert_raises
from clustercorr.feature import ClusterFeature


def test_cluster_feature():

    c = ClusterFeature('chr1', 1, 2, range(5))
    yield check_equal, c.group, 'chr1'
    yield check_equal, c.start, 1
    yield check_equal, c.end, 2
    yield check_equal, c.values, range(5)

def test_cluster_feature_dist():
    c1 = ClusterFeature('chr1', 1, 2, range(5))
    c2 = ClusterFeature('chr1', 1, 2, range(5))

    assert c1.distance(c2) == 0
    c2.start, c2.end = 3, 4
    assert c1.distance(c2) == 1
    c2.start, c2.end = 0, 1000
    assert c1.distance(c2) == 0
    assert c2.distance(c1) == 0


    c2.group = "chr2"
    assert c1.distance(c2) > 1e8


def test_cluster_corr():
    c1 = ClusterFeature('chr1', 1, 2, range(5))
    assert c1.is_correlated(c1)
    c2 = ClusterFeature('chr1', 1, 2, [-x for x in range(5)])
    assert not c1.is_correlated(c2)


def check_equal(a, b, msg=None):
    assert a == b, msg
