
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

options(showWarnCalls=FALSE, warn=-1)

stouffer_liptak = function(pvalues, sigma, lower.tail=TRUE){
    qvalues = qnorm(pvalues, mean=0, sd=1, lower.tail=lower.tail)
    C = chol(sigma)
    Cm1 = solve(C) # C^-1
    qvalues = Cm1 %*% qvalues # Qstar = C^-1 * Q
    Cp = sum(qvalues) / sqrt(length(qvalues))
    pstar = pnorm(Cp, mean=0, sd=1, lower.tail=lower.tail)
    return(list(C=Cp, p=pstar))
}

stouffer_liptak.run = function(covs, formula){
  meth = dcast(covs, CpG ~ id, value.var="methylation")
  rownames(meth) = meth$CpG
  meth = as.matrix(meth[,2:ncol(meth)], nrow=nrow(meth))
  covs = covs[covs$CpG == unique(covs$CpG)[1],]
  mod = model.matrix(formula, covs) 
  covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]

  fit = eBayes(lmFit(meth, mod))
  beta.orig = coefficients(fit)[,covariate]
  pvals = topTable(fit, coef=covariate)[,"P.Value"]
  beta.ave = sum(beta.orig) / length(beta.orig)

  #sigma = cor(t(meth), method="spearman")
  #s = stouffer_liptak(pvals, sigma)
  sigma = cor(t(meth), method="pearson")
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

bumping.run = function(covs, formula, n_sims=100){
    suppressPackageStartupMessages(library('parallel', quietly=TRUE))
    meth = dcast(covs, CpG ~ id, value.var="methylation")
    rownames(meth) = meth$CpG
    meth = as.matrix(meth[,2:ncol(meth)], nrow=nrow(meth))
    acovs = covs[covs$CpG == unique(covs$CpG)[1],]
  
    mod = model.matrix(formula, acovs)
    covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]
    mod0 = mod[,!colnames(mod) == covariate, drop=FALSE]
  
    sim_beta_sums = permute.residuals(meth, mod, mod0, iterations=n_sims)
    stopifnot(length(sim_beta_sums) == n_sims)
    fit = lmFit(meth, mod)
    w = fit$sigma
    icoef = coefficients(fit)[,covariate]
    beta_sum = sum.lowess(icoef, w)

    raw_beta_sum = sum(coefficients(fit)[,covariate])
    ngt = sum(abs(sim_beta_sums) >= abs(beta_sum))
    if(ngt < 4 & n_sims == 100) return(bumping.run(covs, formula, 2000))
    if(ngt < 10 & n_sims == 2000) return(bumping.run(covs, formula, 5000))
    if(ngt < 10 & n_sims == 5000) return(bumping.run(covs, formula, 15000))
    pval = (1 + ngt) / (1 + n_sims)
    return(c(covariate, pval, raw_beta_sum / nrow(meth)))
}

#covs = read.csv(file='t.txt')
#bumping.run(covs, methylation ~ asthma)

gee.run = function(covs, formula, cluster_col="CpG", corstr="ex"){
    suppressPackageStartupMessages(library('geepack', quietly=TRUE))
    stopifnot(!is.null(cluster_col))
    covs$clustervar = covs[,cluster_col]
    s = summary(geeglm(formula, id=clustervar, data=covs, corstr=corstr))

    mm = model.matrix(formula, covs)
    covariate = colnames(mm)[1 + as.integer(colnames(mm)[1] == "(Intercept)")]
    row = s$coefficients[covariate,]
    stat = row[['Wald']]
    return(c(covariate, row[[ 'Pr(>|W|)']], row[['Estimate']]))
}

mixed_model.run = function(covs, formula){
    suppressPackageStartupMessages(library('lme4', quietly=TRUE))
    suppressPackageStartupMessages(library('multcomp', quietly=TRUE))
    m = lmer(formula, covs)
    # take the first column unless it is intercept
    covariate = names(fixef(m))[1 + as.integer(names(fixef(m))[1] == "(Intercept)")]
    s = summary(glht(m, paste(covariate, "0", sep=" == ")))
    return(c(covariate, s$test$pvalues[[1]], s$test$coefficients[[1]]))
}

skat.run = function(covs, formula){
  suppressPackageStartupMessages(library('SKAT', quietly=TRUE))
  meth = dcast(covs, CpG ~ id, value.var="methylation")
  rownames(meth) = meth$CpG
  meth = t(as.matrix(meth[,2:ncol(meth)], nrow=nrow(meth)))
  covariate = all.vars(formula)[1]
  #write.csv(covs, file="/tmp/sk.csv")
  covs = covs[covs$CpG == unique(covs$CpG)[1],]

  capture.output(obj <- SKAT_Null_Model(formula, out_type="D", data=covs))
  #sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="davies", r.corr=0.6, kernel="linear")
  sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="optimal.adj", kernel="linear")
  #sink()
  return(c(covariate, sk$p.value, 'nan'))

}

#covs = read.csv('/tmp/sk.csv')
#skat.run(covs, asthma ~ age + gender + race_white + race_hispanic + race_aa)

clust.lm = function(covs, formula, gee.corstr=NULL, gee.clustervar=NULL, limma.block=NULL, bumping=FALSE, liptak=FALSE, skat=FALSE){
    formula = as.formula(formula)
    stopifnot(is.null(gee.corstr) || is.null(limma.block))
    if(bumping){
        return(bumping.run(covs, formula))
    }
    if(skat){
        return(skat.run(covs, formula))
    }
    if(liptak){
        return(stouffer_liptak.run(covs, formula))
    }
    stopifnot("methylation" %in% colnames(covs))
    is.mixed.model = any(grepl("|", attr(terms(formula), 'term.labels'), fixed=TRUE))
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
        # TODO
        stopifnot(limma.block)
        library('reshape2')
        methylation = melt(covs, id.vars=limma.block, measure.vars="methylation")
    }

}


fclust.lm = function(covs, formula, gee.corstr=NULL, ...){
    covs = read.csv(covs)
    return(clust.lm(covs, formula, gee.corstr=gee.corstr, ...))
}
