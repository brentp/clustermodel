"""
script to generate a bunch of bash/bsub commands so we can run all possible
methods on a cluster in order to compare them.
"""
import os

base_model="methylation ~ case"
covs="work/tcga.fake.covs.txt"
meth="work/tcga.fake.meth.txt"
out = "work/out"


group = os.path.splitext(os.path.basename(covs))[0]
sds = ""
extra = "| bsub -J {name} -e logs/{name}.err -o logs/{name}.out -M 2000000"

base_cmd = "echo 'python -m clustercorr \"{model}\" {covs} {meth} {method} {sds} > {out}/{name}.{group}.bed "
base_cmd += "'" + extra

# switch for skat # methylatoin ~ disease + age => disease + age
sk_model = base_model.split("~")[1].strip()
# disease, +, age
sk_model = sk_model.split()
sk_model = sk_model[0] + " ~ " + (" + ".join(t for t in sk_model[1:] if t != "+") or "1")

for name in ("liptak", "bumping"):
    model = base_model
    method = "--" + name
    print base_cmd.format(**locals())

print base_cmd.format(name="skat", method="--skat", model=sk_model,
        group=group, out=out, sds=sds, covs=covs, meth=meth, extra=extra)

for isds in (0, 3, 4, 5):
    sds = "--outlier-sds %i" % isds
    method = ""

    model = base_model + " + (1|CpG)"
    name = "cpg_intercept-sds_%i" % isds
    print base_cmd.format(**locals())

    model = base_model + " + (1|id)"
    name = "id_intercept-sds_%i" % isds
    print base_cmd.format(**locals())

    model = base_model + " + (1|CpG) + (1|id)"
    name = "both_intercept-sds_%i" % isds
    print base_cmd.format(**locals())

    for gee in ('ar,id', 'ex,CpG'):
        name = "gee-" + gee.replace(',', '-') + ("-sds_%i" % isds)
        model = base_model
        method = "--gee-args " + gee
        print base_cmd.format(**locals())

