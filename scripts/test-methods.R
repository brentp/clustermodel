library(clustermodelr)
library(parallel)
library(data.table)
set.seed(42)

# TODO: test large DMR with few samples and small DMR with 10+ each.
n_sims = 101
n_each = 5 # number of samples from each group
n_probes = 2
rho = 0.21
p_cutoff = 0.05 / 100000 * n_probes # we do fewer tests with more probes
effect = 0.03
sd = 0.01

gen_samples = function(rho, n_each, n_probes, means=c(0, 0.0), sds=c(0.01, 0.01)){
    means = rep(means, each=n_each)
    gen.correlated(rho, n_each * 2, n_probes, mean=means, sd=sds[1])
}

# generate from a t-disribution.
gen_samples.t = function(rho, n_each, n_probes, means=c(0, 0.04), sds=c(0.01, 0.01)){
    cases = means[1] + make.correlated(rho, matrix(sds[1] * rt(n_each * n_probes, 3), nrow=n_each))
    controls = means[2] + make.correlated(rho, matrix(sds[2] * rt(n_each * n_probes, 3), nrow=n_each))
    as.matrix(rbind(cases, controls))
}

covs = data.frame(case=c(rep(1, n_each), rep(0, n_each)))
#print(rowMeans(gen_samples(rho, n_each, n_probes, means=c(0.05, 0))))

#stop()


found.true = mclapply(1:n_sims, function(i){
    meth = gen_samples(rho, n_each, n_probes, means=c(-effect/2, effect/2), sds=c(sd, sd))
    covs.long = expand.covs(covs, meth)
    mm = mixed_modelr(covs.long, methylation ~ case + (1|id) + (1|CpG))$p < p_cutoff
    ge = geer(covs.long, methylation ~ case, corstr="ar", idvar="id")$p < p_cutoff
    ex = geer(covs.long, methylation ~ case, corstr="ex", idvar="id")$p < p_cutoff
    excpg = geer(covs.long, methylation ~ case, corstr="ex", idvar="CpG")$p < p_cutoff
    res = list(ge=ge, mm=mm, ex=ex, excpg=excpg)#, lip=lip, li=li)
    
    res$li = stouffer_liptakr(covs, meth, methylation ~ case, cor.method="spearman")$p < 0.01
    #stouffer_liptakr(covs, meth, methylation ~ case)$p < p_cutoff
    res$lip = stouffer_liptakr(covs, meth, methylation ~ case, cor.method="pearson")$p < 0.01
    res$mmid = mixed_modelr(covs.long, methylation ~ case + (1|id))$p < p_cutoff
    res
}, mc.cores=4)
found.true = rbindlist(found.true)

found.false = mclapply(1:n_sims, function(i){
    meth = gen_samples(rho, n_each, n_probes, means=c(0, 0), sds=c(sd, sd))
    covs.long = expand.covs(covs, meth)
    mm = mixed_modelr(covs.long, methylation ~ case + (1|id) + (1|CpG))$p < p_cutoff
    ge = geer(covs.long, methylation ~ case, corstr="ar", idvar="id")$p < p_cutoff
    ex = geer(covs.long, methylation ~ case, corstr="ex", idvar="id")$p < p_cutoff
    excpg = geer(covs.long, methylation ~ case, corstr="ex", idvar="CpG")$p < p_cutoff
    res = list(ge=ge, mm=mm, ex=ex, excpg=excpg)#, lip=lip, li=li)
    
    res$li = stouffer_liptakr(covs, meth, methylation ~ case, cor.method="spearman")$p < 0.01
    res$lip = stouffer_liptakr(covs, meth, methylation ~ case, cor.method="pearson")$p < 0.01
    res$mmid = mixed_modelr(covs.long, methylation ~ case + (1|id))$p < p_cutoff
    res
}, mc.cores=4)
found.false = rbindlist(found.false)



write("method\ttpr\tfpr\ttpr\ttpr/fpr", stdout())
for(k in colnames(found.true)){
    tpr = sum(found.true[[k]]) / n_sims
    fpr = sum(found.false[[k]]) / n_sims
    write(sprintf("%s\t%.4f\t%.4f\t%.4f", k, tpr, fpr, tpr/fpr)
    , stdout())
}    
