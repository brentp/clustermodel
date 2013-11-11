from nose.tools import assert_raises
from clustercorr.clustermodel import clustered_model_frame
import os.path as op

HERE = op.dirname(__file__)

def test_model_frame():

    fname = op.join(HERE, "example-long.csv")

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
    import pandas as pd
    import tempfile

    fname = op.join(HERE, "example-long.csv")
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


