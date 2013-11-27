from clustermodel.clustermodel import clustered_model
from clustermodel.__main__ import fix_name
import os.path as op
import pandas as pd
import tempfile
import numpy as np

# avoid very verbose reprs with nosetests
pd.DataFrame.__repr__ = lambda self: "<DataFrame>"
pd.Series.__repr__ = lambda self: "<Series>"

HERE = op.dirname(__file__)

def test_model():

    meth = pd.read_csv(op.join(HERE, "example-meth.csv"), index_col=0).T
    meth1 = meth.ix[1, :]

    covs = pd.read_table(op.join(HERE, "example-covariates.txt"))

    for kwargs in ({'gee_args': ('ar', 'id')},
                   {'gee_args': ('ex', 'CpG')},
                   {'combine': 'liptak'},
                   {'combine': 'z-score'},
                   {'bumping': True},):

        yield check_clustered_model, covs, meth, "methylation ~ disease", kwargs
        yield check_clustered_model, covs, meth1, "methylation ~ disease", kwargs


    for random_effects in ("(1|CpG)", "(1|id)", "(1|id) + (1|CpG)"):
        yield (check_clustered_model, covs, meth,
                "methylation ~ disease + " + random_effects, {})

        yield (check_clustered_model, covs, meth1,
                "methylation ~ disease + " + random_effects, {})


def check_clustered_model(covs, meth, model, kwargs):

    res = clustered_model(covs, meth, model, **kwargs)

    #{'p': 0.153760092338262, 'model': 'methylation ~ disease', 'covariate':
    #        'diseaseTRUE', 'liptak': True, 'coef': 0.125455808080808}
    for k in 'p model covariate coef'.split():
        assert k in res, res

    assert ("|" in model) == ("|" in res['model'][0]), (model, res['model'][0])


def test_clustered_model():

    # test for 20 samples and 5 CpGs

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
            r = clustered_model(covs, meth, model, X="'%s'" % fh.name)
            yield check_clustered_df, r, model, exp


def check_clustered(r, model):
    assert 'p' in r
    assert 'coef' in r
    assert isinstance(r['coef'][0], float), (r['coef'], r['coef'][0])
    assert r['model'][0] == model, (model, r['model'][0])


def check_clustered_df(df, model, exp):
    for gene, (i, row) in zip(exp.index, df.iterrows()):
        assert row['X'] == fix_name(gene)
