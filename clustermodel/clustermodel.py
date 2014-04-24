import os
import sys
import warnings
import numpy as np
import pandas as pd
from .pyper import R
from .send_bin import send_arrays

import tempfile
r = R(max_len=5e7, return_err=False)
#r('library(clustermodelr)')
r('source("~/src/clustermodelr/R/clustermodelr.R");source("~/src/clustermodelr/R/combine.R")')
#r('source("/usr/local/src/clustermodelr/R/clustermodelr.R");source("/usr/local/src/clustermodelr/R/combine.R")')

def ilogit(v):
    return 1 / (1 + np.exp(-v))

def kwargs_to_str(kwargs):
    def convert(v):
        if v is True: return "TRUE"
        if v is False: return "FALSE"
        if v is None: return "NA"
        if isinstance(v, basestring):
            return "'%s'" % v
        return v
    return ", ".join("%s=%s"  % (k, convert(v))
                            for k, v in kwargs.iteritems())

def rcall(cov, meths, model, X=None, weights=None, kwargs=None,
        bin_fh=tempfile.NamedTemporaryFile(suffix='.cluster.bin'),
        weight_fh=tempfile.NamedTemporaryFile(suffix='.cluster.bin')):
    """
    internal function to call R and return the result
    """
    if kwargs is None: kwargs = {}

    # send the methylation arrays via binary. this is
    # much faster than relying on pyper to send large
    # matrices. send_arrays does a seek(0).
    send_arrays(meths, bin_fh.file)
    r('meths = read.bin("%s")' % bin_fh.name)
    if weights is not None:

        send_arrays(weights, weight_fh.file)
        r('weights = read.bin("%s")' % weight_fh.name)
    else:
        r('weights = NULL')

    if not 'mc.cores' in kwargs:
        from . import CPUS
        kwargs['mc.cores'] = CPUS

    # faster to use csv than to use pyper's conversion
    # TODO: only send this once.
    if not isinstance(cov, str):
        fh = tempfile.NamedTemporaryFile()
        cov.to_csv(fh, index=False)
        fh.flush()
        cov = fh.name

    assert os.path.exists(cov), cov
    r['cov'] = cov
    if X is None:
        kwargs_str = kwargs_to_str(kwargs)
        #print >>sys.stderr, "fclust.lm(cov, meths, '%s', %s)" % (model, kwargs_str)
        r("a <- data.frame(p=NaN, coef=NaN, covariate=NA); a <- mclust.lm('%s', cov, meths, weights=weights, %s)"
                % (model, kwargs_str))
        df = r['a']
        df['model'] = model
        df['p'] = df['p'].astype(float)
    else:
        kwargs_str = kwargs_to_str(kwargs)
        #print >>sys.stderr, "mclust.lm.X('%s', cov, meths, %s, %s)" % (model, X, kwargs_str)
        r("a = data.frame(p=NaN, coef=NaN, covariate=NA); a <- mclust.lm.X('%s', cov, meths, %s, weights=weights, %s)" % (model, X, kwargs_str))
        df = r['a']

    df['coef'] = df['coef'].astype(float)
    # since we're probably operating on logit transformed data
    # we do the inverse logit and subtract 0.5 since ilogit(0) == 0.5
    df['icoef'] = ilogit(df['coef']) - 0.5
    return df

# TODO: weights should be same shape, type as cluster_dfs here.
def clustered_model(cov_df, cluster_dfs, model, X=None, weights=None, gee_args=(),
        combine=False, bumping=False, skat=False, counts=False, outlier_sds=None):
    """
    Given a cluster of (presumably) correlated CpG's. There are a number of
    methods one could employ to determine the association of the methylation
    of those CpG's to a covariate (often disease status).
    Here we implement:

    1. GEE with autoregressive, independent, or exchangeable correlation
       structure.
    2. random interecpt for each CpG and or each individual
    3. Liptak correction of correlated p-values (1 p-value from each probe and
       use observed correlation for correction
    4. "bumping" algorithm that simulates data by shuffling the residuals of
       the reduced model and compares the observed coefficient estimates to the
       simulated estimates from the shuffled data. Uses a lowess smooothing of
       the data. Different from "bumphunting" because "bumphunting" must do this
       genome-wide, whereas here, it is on a per-cluster basis
    5. SKAT. we can use the methylation values to send to SKAT to compare to the
       null model that does not consider methylation.


    Arguments:

        cov_df - a pandas.DataFrame that must have an index of sample ids and
                 all the covariates defined in model except "id" and "CpG"
                 which will be set automatically. pandas DataFrame

        cluster_df - a pandas.DataFrame that is a cluster of probes. Must have
                     and index indicating the CpG (name or site) and columns
                     of sample ids. This function will use samples from the
                     intersection of cluster_df.columns and cov_df.index

        model - model in R syntax with "methylation ~" as the RHS. Other
                allowed covariates are any that appear in cov_df as well as
                "CpG" and "id" which will be set by this function.
                If not using combine or bumping or gee_args this model should
                likely have a random effect. The obvious choices would be:
                (1|CpG) and/or (1|id) to add a random-intercept by CpG site
                and/or by sample.
                The p-value returned will always be for the first covariate
                in the model. See module docstring for examples.

        X - a file with the same samples as cov_df and rows of expression
            data. If present, each DMR will be tested against each row in
            this file--this is computationally intensive!! Or a data.frame
            or matrix already defined in R.

        gee_args - a 2-tuple of arguments to R's geepack::geeglm().
                   1) the corstr (one of "ex", "in", "ar")
                   2) the cluster variable. This will likely be "id" if
                      corstr="ar" otherwise it will likely be "CpG"
                   sorting is handled internally so that these options will
                   work as expected.
                   So common invocations would be ('ex', 'CpG') or ('ar', 'id')

        combine - either 'liptak' or 'z-score' method for combining the
                  p-values from modelling each CpG independently.

        bumping - if set to True, use a modified bump-hunting algorithm to test
                  the sum of the observed coefficients (1 for each CpG) against
                  the sums of coefficients derived by fitting the model to data
                  generated by repeatedly shuffling the residuals of the reduced
                  model

        skat - if set to True, use skat to test if modelling the CpG
               methylation better describes the dependent variable.
    """

    cov_df['id'] = np.arange(cov_df.shape[0]).astype(int)
    cov = cov_df
    meths = cluster_dfs if not isinstance(cluster_dfs, (pd.DataFrame,
                                                        pd.Series)) \
                        else [cluster_dfs]
    if weights is not None:
        weights = weights if not isinstance(weights, (pd.DataFrame, pd.Series)) \
                      else [weights]

    if outlier_sds > 0:
        [set_outlier_nan(cluster_df, outlier_sds) for cluster_df in meths]

    if "|" in model:
        assert not any((skat, combine, bumping, gee_args))
        return rcall(cov, meths, model, X, weights=weights, kwargs=dict(counts=counts))

    if skat:
        return rcall(cov, meths, model, X, weights=weights, kwargs=dict(skat=True))
    elif combine:
        return rcall(cov, meths, model, X, weights=weights, kwargs=dict(combine=combine))
    elif bumping:
        return rcall(cov, meths, model, X, weights=weights, kwargs=dict(bumping=True))
    elif gee_args:
        corr, col = gee_args
        assert corr[:2] in ('ex', 'ar', 'in', 'un')
        return rcall(cov, meths, model, X, weights=weights,
                kwargs={"gee.corstr": corr, "gee.idvar": col, "counts": counts})
    else:
        raise Exception('must specify one of skat/combine/bumping/gee_args'
                        ' or specify a mixed-effect model in lme4 syntax')

def set_outlier_nan(cluster_df, n_sds):
    """
    take cluster dataframe and set to nan
    any values where that are > n_sds standard-deviations away
    from the mean for that probe
    """
    #imean, isd = cluster_df.mean(axis=1), cluster_df.std(axis=1,
    #        skipna=True)

    for probe in cluster_df.index:
        row = cluster_df.ix[probe, :]
        rown = row[~np.isnan(row)]
        m, s = rown.mean(), rown.std()
        rng = (m - (n_sds * s)), (m + (n_sds * s))
        row[((row < rng[0]) | (row > rng[1]))] = np.nan
