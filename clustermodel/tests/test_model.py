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

def test_weights():

    meth = pd.read_csv(op.join(HERE, "example-meth.csv"), index_col=0).T
    meth1 = meth.ix[1, :]
    weights = pd.DataFrame(meth.copy())
    weights.ix[:, :] = 1.0
    weights1 = weights.ix[1, :]

    weights_diff = weights.copy()
    weights_diff.ix[:, :36] = 10.0


    covs = pd.read_table(op.join(HERE, "example-covariates.txt"))

    for kwargs in ({'gee_args': ('ar', 'id')},
            #       {'gee_args': ('ex', 'CpG')},
                   {'combine': 'liptak'},
                   {'combine': 'z-score'},
                   {'bumping': True},):

        yield check_weights1, covs, meth, weights, "methylation ~ disease", kwargs
        yield check_weights1, covs, meth1, weights1, "methylation ~ disease", kwargs

        yield check_weights_diff, covs, meth, weights_diff, "methylation ~ disease", kwargs

    for random_effects in ("(1|CpG)", "(1|id)", "(1|id) + (1|CpG)"):
        yield (check_weights1, covs, meth, weights,
                "methylation ~ disease + " + random_effects, {})

        yield (check_weights1, covs, meth1, weights1,
                "methylation ~ disease + " + random_effects, {})

        yield (check_weights_diff, covs, meth, weights_diff,
                "methylation ~ disease + " + random_effects, {})

def check_weights_diff(covs, meth, weights, model, kwargs):
    res = clustered_model(covs, meth, model, **kwargs)
    resw = clustered_model(covs, meth, model, weights=weights, **kwargs)

    for k in 'p coef'.split():
        for i in range(len(res[k])):
            assert kwargs.get('bumping') or res[k][i] != resw[k][i], (k, i, res[k][i], resw[k][i])

    for k in 'covariate model'.split():
        for i in range(len(res[k])):
            assert res[k][i] == resw[k][i]

def check_weights1(covs, meth, weights, model, kwargs):

    res = clustered_model(covs, meth, model, **kwargs)
    resw = clustered_model(covs, meth, model, weights=weights, **kwargs)

    for k in 'p model covariate coef'.split():
        assert k in res, res
        assert k in resw, resw
        if not 'bumping' in kwargs or k in ('model', 'covariate'):

            for i in range(len(res[k])):
                val, valw = res[k][i], resw[k][i]
                eq = abs(val - valw) < 1e-4 if isinstance(val, float) \
                                            else val == valw
                assert eq or (np.isnan(res[k][i]) and np.isnan(resw[k][i])),\
                         (res[k][i], resw[k][i], k, i)

    assert ("|" in model) == ("|" in res['model'][0]), (model, res['model'][0])

def check_clustered_model(covs, meth, model, kwargs):

    res = clustered_model(covs, meth, model, **kwargs)

    #{'p': 0.153760092338262, 'model': 'methylation ~ disease', 'covariate':
    #        'diseaseTRUE', 'liptak': True, 'coef': 0.125455808080808}
    for k in 'p model covariate coef'.split():
        assert k in res, res

    assert ("|" in model) == ("|" in res['model'][0]), (model, res['model'][0])

def _make_data(ddiff=0):
    covs = pd.DataFrame({'age': range(1, 21), 'sex': ['M'] * 10 + ['F'] * 10,
        'disease': [True] * 5 + [False] * 5 + [True] * 5 + [False] * 5})
    covs.index = ['sample_%i' % i for i in range(1, 21)]

    meth = pd.DataFrame(np.abs(np.random.randn(5, 20)), index = ["chr1:%i" % (100 * i)
        for i in range(1, 6)])
    meth.columns = list(covs.index)
    meth.ix[:, covs.index[covs.disease]] += ddiff

    return covs, meth

def test_clustered_model():

    # test for 20 samples and 5 CpGs
    covs, meth = _make_data()

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

def test_weighted_model():

    covs, meth = _make_data(1.0)
    weights = meth.copy()
    weights.ix[:, :] = 1
    dis = covs.index[covs.disease == True]
    weights.ix[:, dis] = 0.25

    # so now we have a huge difference, but we down-weight
    # all of the cases, so we should always have a larger p-value
    # with weights
    model = "methylation ~ disease"
    mix_model = "methylation ~ disease + (1|CpG)"
    yield check_weight_m, covs, meth, weights, model, {'combine': 'liptak'}
    yield check_weight_m, covs, meth, weights, model, {'combine': 'z-score'}
    #yield check_weight_m, covs, meth, weights, model, {'gee_args': ('ar', 'id')}
    #yield check_weight_m, covs, meth, weights, model, {'gee_args': ('ex', 'CpG')}
    #yield check_weight_m, covs, meth, weights, model, {'bumping': True}
    #yield check_weight_m, covs, meth, weights, mix_model, {}


    meth.ix[4, 1] = np.nan
    yield check_weight_m, covs, meth, weights, model, {'combine': 'liptak'}
    yield check_weight_m, covs, meth, weights, model, {'combine': 'z-score'}
    #yield check_weight_m, covs, meth, weights, model, {'bumping': True}


def check_weight_m(covs, meth, weights, model, kwargs):

    res = clustered_model(covs, meth, model, **kwargs)
    resw = clustered_model(covs, meth, model, weights=weights, **kwargs)
    import sys
    print sys.stderr, res
    print sys.stderr, resw

    for i in range(len(res)):
        assert resw['p'][i] > res['p'][i], (
                'p', i, resw['p'][i], res['p'][i])
        assert resw['coef'][i] <= res['coef'][i], (
                'coef', i, resw['coef'][i], res['coef'][i])






def check_clustered(r, model):
    assert 'p' in r
    assert 'coef' in r
    assert isinstance(r['coef'][0], float), (r['coef'], r['coef'][0])
    assert r['model'][0] == model, (model, r['model'][0])


def check_clustered_df(df, model, exp):
    for gene, (i, row) in zip(exp.index, df.iterrows()):
        assert row['X'] == fix_name(gene)
