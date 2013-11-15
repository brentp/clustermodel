
# comparisons
# 1. gee with:
#    a. ar1
#    b. exchangeable
# 2. mixed-effects with:
#    a. random-intercept for CpG
#    b. random-intercept by sample_id
#    c. random-intercept by CpG and sample_id
#    d. random-intercept by CpG and sample_id and random slope by CpG
# 3. limma with: (TODO)
#    a. correlation within CpG
#    b. correlation with sample fixed-effect and for CpG
# 4. stouffer-liptak:
#    a. probe-level p-values combined using observed correlation.
# 5. bumphunter
#    a. genome-wide using R package
#    b. local using sum of Betas compared to betas of shuffled residuals

suppressPackageStartupMessages(library('limma', quietly=TRUE))
suppressPackageStartupMessages(library('reshape2', quietly=TRUE))

#options(showWarnCalls=FALSE, warn=-1)

stouffer_liptak = function(pvalues, sigma, lower.tail=TRUE){
    qvalues = qnorm(pvalues, mean=0, sd=1, lower.tail=lower.tail)
    C = chol(sigma)
    Cm1 = solve(C) # C^-1
    qvalues = Cm1 %*% qvalues # Qstar = C^-1 * Q
    Cp = sum(qvalues) / sqrt(length(qvalues))
    pstar = pnorm(Cp, mean=0, sd=1, lower.tail=lower.tail)
    return(list(C=Cp, p=pstar))
}

stouffer_liptak.run = function(covs, meth, formula){
    # TODO: what if missing data in covariates.
    covs$methylation = 1 # 
    sigma = cor(meth, method="spearman")
    meth = t(meth)

    mod = model.matrix(formula, covs)
    covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]

    fit = eBayes(lmFit(meth, mod))
    beta.orig = coefficients(fit)[,covariate]
    pvals = topTable(fit, coef=covariate)[,"P.Value"]
    beta.ave = sum(beta.orig) / length(beta.orig)

    p = stouffer_liptak(pvals, sigma)
    return(c(covariate, p$p, beta.ave))
}

# for bumping
permute.residuals = function(mat, mod, mod0, iterations=100, p_samples=1, mc.cores=10){
    stopifnot(nrow(mod) == ncol(mat))

    reduced_lm = lmFit(mat, mod0)
    reduced_residuals = residuals(reduced_lm, mat)
    reduced_fitted = fitted(reduced_lm)

    fit = lmFit(mat, mod)

    coef.name = setdiff(colnames(mod), colnames(mod0))
    beta.orig = coefficients(fit)[,coef.name]

    rm(reduced_lm, fit); gc()
    nc = ncol(reduced_residuals)

    beta.list = mclapply(1:iterations, function(ix){
        mat_sim = reduced_fitted + reduced_residuals[,sample(1:nc)]
        ifit = lmFit(mat_sim, mod)
        icoef = coefficients(ifit)[,coef.name]
        w = ifit$sigma
        # get names as integer positions:
        names(icoef) = gsub(".*?(\\d+)$", "\\1", names(icoef), perl=T)
        names(w) = gsub(".*?(\\d+)$", "\\1", names(w), perl=T)
        sum.lowess(icoef, w)
    }, mc.cores=mc.cores)
        
    beta.sum = rep(0, n=iterations)
    for(i in 1:iterations){
        beta.sum[i] = beta.list[[i]]
    }
    beta.sum
}

sum.lowess = function(icoefs, weights, span=0.2){
    if(length(icoefs) < 3){ return(sum(icoefs)) }
    res = try(limma::loessFit(icoefs, as.integer(names(icoefs)),
                              span=span, weights=weights), silent=TRUE)
    if(class(res) == "try-error") return(sum(icoefs))
    return(sum(res$fitted))
}

bumping.run = function(covs, meth, formula, n_sims=100){
    suppressPackageStartupMessages(library('parallel', quietly=TRUE))
    covs$methylation = 1 # for formula => model.matrix

    mod = model.matrix(formula, covs)
    covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]
    mod0 = mod[,!colnames(mod) == covariate, drop=FALSE]

    sim_beta_sums = permute.residuals(meth, mod, mod0, iterations=n_sims)
    stopifnot(length(sim_beta_sums) == n_sims)

    fit = lmFit(meth, mod)
    w = fit$sigma

    icoef = coefficients(fit)[,covariate]
    # get names as integer positions:
    names(icoef) = gsub(".*?(\\d+)$", "\\1", names(icoef), perl=T)
    names(w) = gsub(".*?(\\d+)$", "\\1", names(w), perl=T)
    beta_sum = sum.lowess(icoef, w)

    raw_beta_sum = sum(coefficients(fit)[,covariate])
    ngt = sum(abs(sim_beta_sums) >= abs(beta_sum))
    if(ngt < 4 & n_sims == 100) return(bumping.run(covs, meth, formula, 2000))
    if(ngt < 10 & n_sims == 2000) return(bumping.run(covs, meth, formula, 5000))
    if(ngt < 10 & n_sims == 5000) return(bumping.run(covs, meth, formula, 15000))
    pval = (1 + ngt) / (1 + n_sims)
    return(c(covariate, pval, raw_beta_sum / nrow(meth)))
}

#covs = read.csv(file='t.txt')
#bumping.run(covs, methylation ~ asthma)

gee.run = function(covs, formula, cluster_col="CpG", corstr="ex"){
    # assume it's already sorted by CpG, then by id.
    if(cluster_col != "CpG" && corstr == "ar"){
        covs = covs[order(covs[,cluster_col], covs$CpG),]
    }
    suppressPackageStartupMessages(library('geepack', quietly=TRUE))
    stopifnot(!is.null(cluster_col))
    covs$clustervar = covs[,cluster_col]
    # can't do logistc with cluster_col of id, gives bad results for some reason
    s = summary(geeglm(formula, id=clustervar, data=covs, corstr=corstr))
    mm = model.matrix(formula, covs)
    covariate = colnames(mm)[1 + as.integer(colnames(mm)[1] == "(Intercept)")]
    row = s$coefficients[covariate,]
    stat = row[['Wald']]
    return(c(covariate, row[[ 'Pr(>|W|)']], row[['Estimate']]))
}

#gee.run(read.csv('tt.csv'), methylation ~ disease, "id", "ex")

mixed_model.run = function(covs, formula){
    suppressPackageStartupMessages(library('lme4', quietly=TRUE))
    suppressPackageStartupMessages(library('multcomp', quietly=TRUE))
    # automatically do logit regression.
    m = lmer(formula, covs)
    # take the first column unless it is intercept
    covariate = names(fixef(m))[1 + as.integer(names(fixef(m))[1] == "(Intercept)")]
    s = summary(glht(m, paste(covariate, "0", sep=" == ")))
    return(c(covariate, s$test$pvalues[[1]], s$test$coefficients[[1]]))
}


skat.run = function(covs, meth, formula, r.corr=c(0.00, 0.015, 0.06, 0.15)){
    suppressPackageStartupMessages(library('SKAT', quietly=TRUE))
    covariate = all.vars(formula)[1]

    capture.output(obj <- SKAT_Null_Model(formula, out_type="D", data=covs))
    #sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="davies", r.corr=0.6, kernel="linear")
    sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="optimal.adj", kernel="linear",
            r.corr=r.corr)
    #sk <- SKAT(meth, obj, is_check_genotype=TRUE, method="optimal.adj", kernel="linear.weighted", weights.beta=c(1, 10))
    #sk <- SKAT(meth, obj, is_check_genotype=TRUE, method="optimal.adj", kernel="linear")
    #sink()
    return(c(covariate, sk$p.value, 'nan'))
}

#covs = read.csv('/tmp/sk.csv')
#skat.run(covs, asthma ~ age + gender + race_white + race_hispanic + race_aa)

long.covs = function(covs, meth){
    n_samples = nrow(covs)
    stopifnot(nrow(meth) == n_samples)
    # e.g. meth is 68 patients * 4 CpGs
    #      covs is 68 patients * 5 covariates
    # need to replicated covs 4 times (1 per CpG)
    covs = covs[rep(1:nrow(covs), ncol(meth)),]
    cpgs = gsub("CpG__", '', colnames(meth), fixed=TRUE)
    dim(meth) = NULL
    covs$methylation = meth
    covs$CpG = rep(cpgs, each=n_samples) # 1 1 1, 2 2 2, etc since CpG's are grouped.
    covs
}


clust.lm = function(covs, formula, meth=NULL,
                    gee.corstr=NULL, gee.clustervar=NULL,
                    limma.block=NULL, bumping=FALSE, liptak=FALSE, skat=FALSE){

    formula = as.formula(formula)
    stopifnot(is.null(gee.corstr) || is.null(limma.block))

    # we assume there is one extra column for each CpG
    if(is.null(meth)){
        meth = as.matrix(covs[,grep("CpG__", colnames(covs), fixed=TRUE)])
    }
    rownames(meth) = rownames(covs)

    if(bumping){ # wide
        return(bumping.run(covs, t(meth), formula))
    }
    if(skat){ # wide
        return(skat.run(covs, meth, formula))
    }
    if(liptak){ # wide
        return(stouffer_liptak.run(covs, meth, formula))
    }

    ###########################################
    # GEE and mixed models require long format.
    ###########################################
    covs = long.covs(covs[,grep("CpG__", colnames(covs), fixed=TRUE, invert=TRUE)], meth)
    is.mixed.model = any(grepl("|", attr(terms(formula), 'term.labels'), fixed=TRUE))
    #write.csv(covs, 'long.covs.t.csv')
    # mixed-model
    if (is.null(gee.corstr) && is.null(limma.block)){
        stopifnot(is.mixed.model)
        return(mixed_model.run(covs, formula))
    # GEE
    } else if (!is.null(gee.corstr)){
        stopifnot(!is.null(gee.clustervar))
        return(gee.run(covs, formula, cluster_col=gee.clustervar, corstr=gee.corstr))
    # limma
    } else {
        # TODO this goes in the matrix section above and uses
        # duplicateCorrelation
        stop()
    }

}



fclust.lm = function(covs, formula, gee.corstr=NULL, ..., mc.cores=4){

    if(is.character(covs)) covs = read.csv(covs)
    if(!"cluster_set" %in% colnames(covs)){
        res = (clust.lm(covs, formula, gee.corstr=gee.corstr, ...))
        return(data.frame(list(covariate=res[1], p=res[2], coef=res[3])))
    }
    suppressPackageStartupMessages(library('data.table', quietly=TRUE))
    suppressPackageStartupMessages(library('parallel', quietly=TRUE))

    cluster_ids = covs$cluster_set[!duplicated(covs$cluster_set)]
    results = mclapply(cluster_ids, function(cs){
        res = clust.lm(covs[covs$cluster_set == cs,], formula, gee.corstr=gee.corstr, ...)
        list(covariate=res[1], p=res[2], coef=res[3], cluster_id=cs)
    }, mc.cores=mc.cores)
    results = rbindlist(results)
    rownames(results) = cluster_ids
    results
}

if(FALSE){
    covs = read.csv('covs.wide.txt')
    #print(clust.lm(covs, methylation ~ gene.E + disease, bumping=TRUE))
    print(fclust.lm(covs, methylation ~ gene.E + disease + (1|id) + (1|CpG)))
    print(clust.lm(covs, methylation ~ gene.E + disease + (1|id)))
    #print(clust.lm(covs, methylation ~ gene.E, gee.clustervar="id", gee.corstr="ex"))
    print(clust.lm(covs, methylation ~ gene.E, gee.clustervar="id", gee.corstr="ex"))
    print(clust.lm(covs, methylation ~ gene.E, gee.clustervar="id", gee.corstr="ar"))
    print('liptak')
    print(clust.lm(covs, methylation ~ gene.E, liptak=TRUE))
    print('bumping')
    print(clust.lm(covs, methylation ~ gene.E, bumping=TRUE))
    #print(clust.lm(covs, disease ~ gene.E, skat=TRUE))
    #print(clust.lm(covs, methylation ~ gene.E, gee.clustervar="CpG", gee.corstr="un"))
}


fclust.lm.X = function(covs, formula, X, gee.corstr=NULL, ..., mc.cores=12, testing=FALSE){
    library(parallel)
    library(data.table)
    formula = as.formula(formula)
    if(is.character(covs)) covs = read.csv(covs)

    X = read.delim(gzfile(X), row.names=1)
    mc.cores = min(mc.cores, ncol(X))

    if(testing) X = X[400:408,]

    stopifnot(nrow(covs) %% ncol(X) == 0)
    n_each = nrow(covs) / ncol(X)
    # need this for when X_locs is not specified since we never readi
    # in the array in python
    rownames(X) = gsub("-|:| ", ".", as.character(rownames(X)))
    rnames = rownames(X)

    # get a + b + c from y ~ a + b + x
    rhs = as.character(formula)[length(as.character(formula))]
    lhs = as.character(formula)[2]
    irows = 1:nrow(X)
    stopifnot(n_each >= 1)

    results = mclapply(irows, function(irow){
        row = rep(t(X[irow,]), each=n_each)
        covs2 = covs # make a copy so we dont end up with huge covs
        # add the expression column to the dataframe.
        covs2[,rnames[irow]] = row
        sformula = sprintf("%s ~ %s + %s", lhs, rnames[irow], rhs)
        # this doesnt work!!
        #sformula = sprintf("%s ~ %s + %s", rnames[irow], lhs, rhs)
        res = try(clust.lm(covs2, as.formula(sformula),
                           gee.corstr=gee.corstr, ...), silent=FALSE)

        #res.df[irow,1:4] = c(res[1], res[2], res[3], rnames[irow])
        if(!inherits(res, "try-error")){
            if(is.na(res[2])){ stop(res) }
            return(list(covariate=res[1], p=res[2], coef=res[3],
                              X=rnames[irow], model=sformula))
        }else{
            return(list(covariate="error", p=NaN, coef=NaN, X=rnames[irow],
                              model=sformula))
        }

    }, mc.cores=mc.cores)
    results = rbindlist(results)
    rownames(results) = results$X
    results
}

cprint = function(...) write(..., stdout())

test_X = function(){
    covs = "clustercorr/tests/example-wide.csv"
    covs = "covs.wide.multi.csv"
    X = 'clustercorr/tests/example-expression.txt.gz'

    cprint("\nmixed-effects model")
    formula = methylation ~ disease + (1|id) + (1|CpG)
    df = fclust.lm.X(covs, formula, X, testing=TRUE)
    write.csv(df, 'Xout.csv')
    stop()
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nGEE")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = fclust.lm.X(covs, formula, X, testing=TRUE, gee.corstr="ar", gee.clustervar="id")
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nbumping")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = fclust.lm.X(covs, formula, X, testing=TRUE, bumping=TRUE)
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nliptak")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = fclust.lm.X(covs, formula, X, testing=TRUE, bumping=TRUE)
    print(head(df[order(as.numeric(df$p)),], n=5))
    # show that we get the same result (about with the linear model)
    # pvalue is  2.85844757130782e-06 for the clustered approach and
    # 7.88e-07 for looking at a single probe with a linear model in
    # the region. coefficients vary by ~ 0.001.
    probe = "A_33_P3403576"
    X = read.delim(gzfile(X), row.names=1, nrows=408)
    covs = read.csv(covs)
    #covs = covs[covs$CpG == covs$CpG[1],]
    covs$X = t(X[probe,])
    for(cidx in grep("CpG__", colnames(covs), fixed=TRUE)){
        covs$methylation = covs[,cidx]
        cprint(paste0("\n", probe))
        print(summary(lm(methylation ~ X + disease, covs))$coefficients)
    }

}
test_X()
