from nose.tools import assert_raises
from clustercorr.clustermodel import clustered_model_frame, clustered_model
from clustercorr.__main__ import fix_name
import os.path as op
import pandas as pd
import tempfile

HERE = op.dirname(__file__)

def test_model_frame():

    fname = op.join(HERE, "example-wide.csv")

    for kwargs in ({'gee_args': ('ar', 'id')},
                   {'gee_args': ('ex', 'CpG')},
                   {'liptak': True},
                   {'bumping': True},
                   ):

        yield check_clustered_model_frame, fname, "methylation ~ disease", kwargs

    for random_effects in ("(1|CpG)", "(1|id)", "(1|id) + (1|CpG)"):
        yield (check_clustered_model_frame, fname,
                "methylation ~ disease + " + random_effects, {})

def test_model_frame_X():

    fname = op.join(HERE, "example-wide.csv")
    X = op.join(HERE, "example-expression.txt.gz")

    # get small subset for test
    df = pd.read_table('clustercorr/tests/example-expression.txt.gz',
             compression='gzip', index_col=0)

    with tempfile.NamedTemporaryFile() as fh:
        sub = df.ix[:10, :]
        sub.to_csv(fh, index=True, index_label="probe",
                            sep="\t", quote=False)

        fh.flush()
        res = clustered_model_frame(fname,
                "methylation ~ disease + (1|CpG)", X=fh.name)

        assert (res.index == sub.index).all()
        for col in "covariate p coef X model".split():
            col in res.columns, res.columns



def check_clustered_model_frame(fname, model, kwargs):

    res = clustered_model_frame(fname, model, **kwargs)

    #{'p': 0.153760092338262, 'model': 'methylation ~ disease', 'covariate':
    #        'diseaseTRUE', 'liptak': True, 'coef': 0.125455808080808}
    for k in 'p model covariate coef'.split():
        assert k in res, res

    if 'liptak' in kwargs:
        assert res['liptak']
    elif 'gee_args' in kwargs:
        assert 'gee.corstr' in res, res
        assert not 'liptak' in res, res
        assert not 'bumping' in res, res

    elif 'bumping' in res:
        assert res['bumping']
        assert not 'liptak' in res, res

    else: # mixed model
        assert "|" in res['model']
        assert not 'liptak' in res, res
        assert not 'bumping' in res, res


def test_clustered_model():

    # test for 20 samples and 5 CpGs
    import numpy as np

    covs = pd.DataFrame({'age': range(1, 21), 'sex': ['M'] * 10 + ['F'] * 10,
        'disease': [True] * 5 + [False] * 5 + [True] * 5 + [False] * 5})
    covs.index = ['sample_%i' % i for i in range(1, 21)]

    meth = pd.DataFrame(np.random.randn(5, 20), index = ["chr1:%i" % (100 * i)
        for i in range(1, 6)])
    meth.columns = list(covs.index)

    model = "methylation ~ disease + (1|id)"

    r = clustered_model(covs, meth, model)
    yield check_clustered, r, model

    np.random.seed(42)
    exp = meth.copy() * 1.15 + np.random.random(meth.shape)
    for bad_name in ("", "-", " "):
        with tempfile.NamedTemporaryFile(delete=True) as fh:
            exp.index = ['gene' + bad_name + l for l in 'ABCDE']
            exp.to_csv(fh.name, sep="\t", quote=False, index=True,
                    index_label="probe")
            fh.flush()
            r = clustered_model(covs, meth, model, X=fh.name)
            yield check_clustered_df, r, model, exp


def check_clustered(r, model):
    assert 'p' in r
    assert 'coef' in r
    assert isinstance(r['coef'], float)
    assert r['model'] == model


def check_clustered_df(df, model, exp):
    for gene, (i, row) in zip(exp.index, df.iterrows()):
        assert row['X'] == fix_name(gene)
