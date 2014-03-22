library(clustermodelr)
library(parallel)
library(data.table)
set.seed(42)

# TODO: test large DMR with few samples and small DMR with 10+ each.
n_sims = 1530
n_each = 5 # number of samples from each group
n_probes = 2
rho = 0.21
p_cutoff = 0.05 / 5000 * n_probes # we do fewer tests with more probes
effect = 0.015
sd = 0.015

gen_samples = function(rho, n_each, n_probes, means=c(0, 0.0), sds=c(0.01, 0.01)){
    means = rep(means, each=n_each)
    gen.correlated(rho, n_each * 2, n_probes, mean=means, sd=sds[1])
}

# generate from a t-disribution.
gen_samples.t = function(rho, n_each, n_probes, means=c(0, 0.04), sds=c(0.01, 0.01), df=5){
    cases = means[1] + make.correlated(rho, matrix(sds[1] * rt(n_each * n_probes, df), nrow=n_each))
    controls = means[2] + make.correlated(rho, matrix(sds[2] * rt(n_each * n_probes, df), nrow=n_each))
    as.matrix(rbind(cases, controls))
}

covs = data.frame(case=c(rep(1, n_each), rep(0, n_each)))
#print(rowMeans(gen_samples(rho, n_each, n_probes, means=c(0.05, 0))))

#stop()
name.it = function(name, i, res){
    res$name = name
    res$i = i
    res
}


found.true = mclapply(1:(n_sims/20), function(i){
    meth = gen_samples(rho, n_each, n_probes, means=c(-effect/2, effect/2), sds=c(sd, sd))
    covs.long = expand.covs(covs, meth)
    res = rbindlist(list(
          name.it("mm", i, mixed_modelr(covs.long, methylation ~ case + (1|id) + (1|CpG))),
          name.it("gear", i, geer(covs.long, methylation ~ case, corstr="ar", idvar="id")),
          name.it("geex", i, geer(covs.long, methylation ~ case, corstr="ex", idvar="id")),
          name.it("gecpg", i, geer(covs.long, methylation ~ case + id, corstr="ex", idvar="CpG")),
          name.it("lis", i, stouffer_liptakr(covs, meth, methylation ~ case, cor.method="spearman")),
          name.it("lip", i, stouffer_liptakr(covs, meth, methylation ~ case, cor.method="pearson"))
    ))
    res
}, mc.cores=4)

found.true = data.frame(rbindlist(found.true))
found.true$effect = effect


found.false = mclapply(1:n_sims, function(i){
    meth = gen_samples(rho, n_each, n_probes, means=c(0, 0), sds=c(sd, sd))
    covs.long = expand.covs(covs, meth)
    res = rbindlist(list(
          name.it("mm", i, mixed_modelr(covs.long, methylation ~ case + (1|id) + (1|CpG))),
          name.it("gear", i, geer(covs.long, methylation ~ case, corstr="ar", idvar="id")),
          name.it("geex", i, geer(covs.long, methylation ~ case, corstr="ex", idvar="id")),
          name.it("gecpg", i, geer(covs.long, methylation ~ case + id, corstr="ex", idvar="CpG")),
          name.it("lis", i, stouffer_liptakr(covs, meth, methylation ~ case, cor.method="spearman")),
          name.it("lip", i, stouffer_liptakr(covs, meth, methylation ~ case, cor.method="pearson"))
    ))
    res
}, mc.cores=4)
found.false = data.frame(rbindlist(found.false))
found.false$effect = 0

found = rbind(found.true, found.false)
write.table(found, row.names=FALSE, quote=FALSE, file="cluster.cmp.txt", sep="\t")


write("method\ttpr\tfpr\tppv\tprecision\tf1.score", stdout())

for(method in unique(found$name)){
    trues = found[(found$name == method) & (found$effect != 0),]
    falses = found[(found$name == method) & (found$effect == 0),]

    tp = sum(trues$p < p_cutoff)
    fn = nrow(trues) - tp

    fp = sum(falses$p < p_cutoff)
    # recall == tpr
    # precision == number of correct results divided by the number of results that should have been returned.
    tpr = tp / nrow(trues)
    fpr = fp / nrow(falses)

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1.score = 2 * (precision * recall) / (precision + recall)

    write(sprintf("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f", method, tpr, fpr, tpr/(tpr + fpr),
                  precision, f1.score)
    , stdout())
}    
