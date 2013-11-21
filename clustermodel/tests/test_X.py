from nose.tools import assert_raises
from clustermodel.__main__ import distX


def test_distX():

    expr = {'chrom': 'chrA', 'start': 2, 'end': 3, 'strand': '+'}

    for strand in ('+', '-'):
        dmr = {'chrom': 'chrA', 'start': 1, 'end': 2}
        expr['strand'] = strand
        distX(dmr, expr)
        assert dmr['distance'] == 0, (dmr, expr)

    for attr in 'start end strand'.split():
        assert dmr['X' + attr] == expr[attr]


def test_distX_gt0():
    dmr = {'chrom': 'chrA', 'start': 1, 'end': 2}
    expr = {'chrom': 'chrA', 'start': 4, 'end': 5, 'strand': '+'}

    for strand in ('+', '-'):
        expr['strand'] = strand
        distX(dmr, expr)

        if strand == "+":
            assert dmr['distance'] == -2
        else:
            assert dmr['distance'] == 2

