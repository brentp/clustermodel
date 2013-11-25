"""
script to generate a bunch of bash/bsub commands so we can run all possible
methods on a cluster in order to compare them.
"""
import os
import sys
import pandas as pd
import numpy as np

np.random.seed(42)

def shuffle_expr(fexpr):
    """
    break the relation between expression and methylation.
    """
    df = pd.read_table(fexpr, index_col=0)
    orig_cols = list(df.columns)
    cols = np.array(df.columns)
    np.random.shuffle(cols)
    df.columns = cols
    new_f = fexpr.replace('.txt', '.shuffled.txt')
    df[orig_cols].to_csv(new_f, sep="\t", index=True,
            quote=False, float_format="%.4f")
    return new_f


base_model = "methylation ~ 1"

covs = sys.argv[1]
meth = sys.argv[2]

orig_expr = "/proj/Schwartz/brentp/2013/tcga-methex/brca-expr.matrix.txt"
fake_expr = shuffle_expr(orig_expr)


expr_locs = "/proj/Schwartz/brentp/2013/tcga-methex/expr-probe-locs.bed"


group = os.path.splitext(os.path.basename(covs))[0]
sds = ""

extra = "| bsub -J {name} -e logs/{name}.err -o logs/{name}.out -M 2000000"

base_cmd = ("echo 'python -m clustermodel \"{model}\" {covs} {meth} {method} {sds}"
           " --X {expr} --X-locs {expr_locs} --X-dist 50000"
           " > {out}/{name}.{group}.pvals.bed "
           "'" + extra)


# switch for skat # methylatoin ~ disease + age => disease + age
sk_model = base_model.split("~")[1].strip()
# disease, +, age
sk_model = sk_model.split()
sk_model = sk_model[0] + " ~ " + (" + ".join(t for t in sk_model[1:] if t != "+") or "1")

for expr in (orig_expr, fake_expr):
    out = sys.argv[3] + ("fake/" if "shuff" in expr else "real/")
    try:
        os.makedirs(out)
    except OSError:
        pass

    for name in ("liptak", "bumping"):
        model = base_model
        method = "--" + name
        print base_cmd.format(**locals())

    print base_cmd.format(name="skat", method="--skat", model=sk_model,
            group=group, out=out, sds=sds, covs=covs, meth=meth, extra=extra,
            expr=expr, expr_locs=expr_locs)

    for isds in (0, 3):
        sds = "--outlier-sds %i" % isds
        method = ""

        model = base_model + " + (1|CpG) + (1|id)"
        name = "both_intercept-sds_%i" % isds
        print base_cmd.format(**locals())

        for gee in ('ar,id', 'ex,id'):
            name = "gee-" + gee.replace(',', '-') + ("-sds_%i" % isds)
            model = base_model
            method = "--gee-args " + gee
            print base_cmd.format(**locals())

