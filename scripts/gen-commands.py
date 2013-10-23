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
extra = "| bsub -J {name} -e logs/{name}.err -o logs/{name}.out -M 2"

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

method = ""
for isds in (0, 3, 4, 5):
    sds = "--outlier-sds %i" % isds

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


"""
echo "commands=(
python -m clustercorr $skat_model $covs $meth --skat         > $out/skat.$(basename $covs .txt).bed
python -m clustercorr \"$base_model\" $covs $meth --liptak   > $out/liptak.$(basename $covs .txt).bed
python -m clustercorr \"$base_model\" $covs $meth --bumping  > $out/bumping.$(basename $covs .txt).bed

python -m clustercorr \"$base_model + (1|CpG)\" $covs $meth  > $out/cpg_intercept.$(basename $covs .txt).bed
python -m clustercorr \"$base_model + (1|id)\" $covs $meth   > $out/id_intercept.$(basename $covs .txt).bed
python -m clustercorr \"$base_model + (1|CpG) + (1|id)\" $covs $meth > $out/cpg_intercept_id_intercept.$(basename $covs .txt).bed

python -m clustercorr \"$base_model\" $covs $meth --gee-args ar,id  > $out/gee-ar_id.$(basename $covs .txt).bed
python -m clustercorr \"$base_model\" $covs $meth --gee-args ex,CpG > $out/gee-ex_CpG.$(basename $covs .txt).bed

"

# outliers!!
for SDS in 3 4 5; do
    echo "

python -m clustercorr \"$base_model\" + (1|CpG) $covs $meth --outlier-sds $SDS > $out/cpg_intercept.$(basename $covs .txt).outlier-sds-$SDS.bed
python -m clustercorr \"$base_model\" + (1|id) $covs $meth --outlier-sds $SDS  > $out/id_intercept.$(basename $covs .txt).outlier-sds-$SDS.bed
python -m clustercorr \"$base_model\" + (1|CpG) + (1|id)\" $covs $meth --outlier-sds $SDS > $out/cpg_intercept_id_intercept.$(basename $covs .txt).outlier-sds-$SDS.bed

python -m clustercorr \"$base_model\" $covs $meth --gee-args ar,id --outlier-sds $SDS > $out/gee-ar_id.$(basename $covs .txt).outlier-sds-$SDS.bed
python -m clustercorr \"$base_model\" $covs $meth --gee-args ex,CpG --outlier-sds $SDS > $out/gee-ex_CpG.$(basename $covs .txt).outlier-sds-$SDS.bed
"
done

echo ")"
"""
