library(clustermodelr)
library(parallel)
library(data.table)
set.seed(42)

# TODO: test large DMR with few samples and small DMR with 10+ each.
n_sims = 201
n_each = 3 # number of samples from each group
n_probes = 15
rho = 0.2
p_cutoff = 0.005
effect = 0.02

gen_samples = function(rho, n_each, n_probes, means=c(0, 0.0), sds=c(0.02, 0.02)){
    cases = gen.correlated(rho, n_each, n_probes, mean=means[1], sd=sds[1])
    controls = gen.correlated(rho, n_each, n_probes, mean=means[2], sd=sds[2])
    as.matrix(rbind(cases, controls))
}

# generate from a t-disribution.
gen_samples = function(rho, n_each, n_probes, means=c(0, 0.04), sds=c(0.02, 0.02)){
    cases = make.correlated(rho, matrix(sds[1] * rt(n_each * n_probes, 3) + means[1], nrow=n_each))
    controls = make.correlated(rho, matrix(sds[2] * rt(n_each * n_probes, 3) + means[2], nrow=n_each))
    as.matrix(rbind(cases, controls))
}

covs = data.frame(case=c(rep(1, n_each), rep(0, n_each)))

#meth = gen_samples(rho, n_each, 200, means=c(0, 0.0), sds=c(0.03, 0.03))
#print(mean(unlist(lapply(1:(ncol(meth) - 1), function(i) cor(meth[,i], meth[,i + 1] + 1)))))



found = mclapply(1:n_sims, function(i){
    meth = gen_samples(rho, n_each, n_probes, means=c(0, effect), sds=c(0.02, 0.02))
    covs.long = expand.covs(covs, meth)
    mm = mixed_modelr(covs.long, methylation ~ case + (1|id) + (1|CpG))$p < p_cutoff
    ge = geer(covs.long, methylation ~ case, corstr="ar", idvar="id")$p < p_cutoff
    ex = geer(covs.long, methylation ~ case, corstr="ex", idvar="id")$p < p_cutoff
    excpg = geer(covs.long, methylation ~ case, corstr="ex", idvar="CpG")$p < p_cutoff
    
    #li = stouffer_liptakr(covs, meth, methylation ~ case)$p < p_cutoff
    #lip = stouffer_liptakr(covs, meth, methylation ~ case, cor.method="pearson")$p < p_cutoff
    list(ge=ge, mm=mm, ex=ex, excpg=excpg)#, lip=lip, li=li)
}, mc.cores=1)

df = rbindlist(found)


for(k in colnames(df)){
    write(sprintf("%s\t%.4f", k, sum(df[[k]]) / n_sims), stdout())
}    
