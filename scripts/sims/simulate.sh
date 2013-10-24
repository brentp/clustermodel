
linkage=complete
meth=work/tcga.matrix.txt
base=work/simulated

for ns in 10 20 40; do
    for rho in 0.3 0.6; do
        for w in 0 0.8; do
            name=ns_${ns}-rho_${rho}-w_${w}

            prefix=$base/ns_${ns}/rho_${rho}/w_${w}
            mkdir -p $prefix

            gen_cmd="python -m clustercorr simulate --n-samples $ns --rho-min $rho -w $w --linkage $linkage $meth $prefix/sim"
            fit_cmd="python scripts/gen-commands.py $prefix/sim.covs.txt $prefix/sim.meth.txt $prefix/ | bash"
            echo "set -xe; $gen_cmd ; $fit_cmd" \
                | bsub -J $name -e logs/$name.err -o logs/$name.out
        done
    done
done
