clustermodel
============

fitting models to clustered correlated data.

The Problem
===========

You have methylation data (or other locally correlated data) with small,
subtle changes that can not be detected at the single-CpG level and you
want to detect regional changes.  Here is an example of the subtle changes
across a region with disease in green vs control in blue (image created
by this software):

![DMR Plot](https://gist.github.com/brentp/204899d3d220ce8a6acd/raw/dd2ef05d22011a77f8020377eede94942ec6e01c/dmr.png)

The Solution(s)
===============

There are a number of methods for **combining** probes. The difficulty is that
we can't use a simple linear model because we have **correlated** data. The most
obvious approaches are:

 1) **bumphuting** from Irizarry's group: find regions where a covariate from a
    model is consistently higher than by chance. In general "chance" is determined
    by shuffling the residuals of a reduced model, and re-calculating simulated
    betas. (they also implement smoothing and a number of nice additional features)

 2) **mixed-effect models**: fit a random intercept per CpG or per sample to
    account for correlated data.

 3) **GEE models**: use generalized estimating equations. This allows us,
    for example to fit an auto-regressive correlation structure, such as we
    might expect to see between adjacent probes. This was implemented nicely
    in `A-clust` (for exchange-able correlation structure only).

 4) **stouffer-liptak** calculate the p-values for adjacent probes
    independently then combine the p-values, taking the correlation
    among the probes into account.

This package attempts to implement (some form of) all of those so that we can
model our data with any one of them with a small change of parameters. To do
so, it first **clusters the data** using `aclust`.

Tutorial
========

For the command-line, we assume the data is set up so that we have a file of
covariates with rows of samples and a file of methylation data that is of shape
**n_probes X n_samples**.

For example, we may wish to use the stouffer-liptak correction on our data in
this repository to find regions where methylation is related to disease. We use the
model in R syntax: *methylation ~ disease* where the first covariate on the RHS is
always our covariate of interest (where we get the p-value).

```Shell
python -m clustermodel 'methylation ~ disease' \
          clustermodel/tests/example-covariates.txt \
          clustermodel/tests/example-methylation.txt.gz \
          --combine liptak \
          --min-cluster-size 3 | head
```

The output from that looks like this:

	#chrom	start	end	coef	p	n_probes	model	method
	chr1	2043760	2043853	0.0560277777778	0.576888128824	3	methylation ~ disease	liptak
	chr1	2058229	2059086	-0.0795509259259	0.155520943237	6	methylation ~ disease	liptak
	chr1	2063798	2064441	-0.042287037037	0.508509656661	6	methylation ~ disease	liptak
	chr1	2066005	2066446	-0.0178055555556	0.102638244865	4	methylation ~ disease	liptak
	chr1	2081983	2082349	-0.0995648148148	0.103685852362	3	methylation ~ disease	liptak
	chr1	2084318	2084595	-0.0962291666667	0.119301491347	4	methylation ~ disease	liptak
	chr1	2100231	2100435	-0.0035	0.222498206941	3	methylation ~ disease	liptak
	chr1	2106221	2106620	-0.0490555555556	0.424566020253	4	methylation ~ disease	liptak
	chr1	2117424	2118131	-0.029962962963	0.421006269346	3	methylation ~ disease	liptak

where the `p` column is for the combined probes. We can repeat the same using a mixed model with
a random intercept by `CpG` (this covariate becomes available automatically if your data
follows our conventions (see example data for details).

```Shell
python -m clustermodel 'methylation ~ disease + (1|CpG)' \
          clustermodel/tests/example-covariates.txt \
          clustermodel/tests/example-methylation.txt.gz \
          --min-cluster-size 3 | head

#chrom  start   end coef    p   n\_probes    model   method
chr1    2043760 2043853 0.0560277777778 0.385222268513  3   methylation ~ disease + (1|CpG) mixed-model
chr1    2058229 2059086 -0.0795509259259    0.018857108935  6   methylation ~ disease + (1|CpG) mixed-model
chr1    2063798 2064441 -0.042287037037 0.157777826794  6   methylation ~ disease + (1|CpG) mixed-model
chr1    2066005 2066446 -0.0178055555556    0.370739004967  4   methylation ~ disease + (1|CpG) mixed-model
chr1    2081983 2082349 -0.0995648148148    0.0460294040919 3   methylation ~ disease + (1|CpG) mixed-model
chr1    2084318 2084595 -0.0962291666667    0.0115980612503 4   methylation ~ disease + (1|CpG) mixed-model
chr1    2100231 2100435 -0.0035 0.8958980765    3   methylation ~ disease + (1|CpG) mixed-model
chr1    2106221 2106620 -0.0490555555556    0.14617070189   4   methylation ~ disease + (1|CpG) mixed-model
chr1    2117424 2118131 -0.029962962963 0.195491513628  3   methylation ~ disease + (1|CpG) mixed-model
```

Note that we change the model to add the random intercept in `lme4` syntax.
For most cases, we have similar p-values (and identical coefficients) as we did
for the liptak method.

clustering
----------

If we run the same model with different clustering parameters, we may get different resuts:

```Shell
python -m clustermodel 'methylation ~ disease + (1|CpG)' \
          clustermodel/tests/example-covariates.txt \
          clustermodel/tests/example-methylation.txt.gz \
          --min-cluster-size 2 \
          --linkage single  \
          --rho-min 0.2 | head

#chrom  start   end coef    p   n\_probes    model   method
chr1    2043438 2043853 0.0254555555556 0.536712997985  5   methylation ~ disease + (1|CpG) mixed-model
chr1    2046505 2046928 0.0907055555556 9.47135696294e-06   5   methylation ~ disease + (1|CpG) mixed-model
chr1    2058229 2059086 -0.0795509259259    0.018857108935  6   methylation ~ disease + (1|CpG) mixed-model
chr1    2060013 2060642 -0.0391574074074    0.165767587931  3   methylation ~ disease + (1|CpG) mixed-model
chr1    2063798 2064441 -0.042287037037 0.157777826794  6   methylation ~ disease + (1|CpG) mixed-model
chr1    2064764 2064819 0.0565277777778 0.107191645281  2   methylation ~ disease + (1|CpG) mixed-model
chr1    2065241 2065263 -0.0539583333333    0.163734202808  2   methylation ~ disease + (1|CpG) mixed-model
chr1    2066005 2066981 -0.0329768518519    0.0612285551598 6   methylation ~ disease + (1|CpG) mixed-model
chr1    2081983 2082522 -0.0603722222222    0.149962950932  5   methylation ~ disease + (1|CpG) mixed-model
```

Existing Regions
================
We may have a list of regions from one study to compare to another study. We
can do this by specifying the --regions arg that gives a BED file of regions
to test. In this case, the clustering is not performed (so clustering
parameters are not used. Probes that fall into each region in the
--regions argument will be assumed as a cluster and the linear model will be run.
Here is an example call:

    python -m clustermodel \
        'methylation ~ disease + (1|CpG)' \
        clustermodel/tests/example-covariates.txt \
        clustermodel/tests/example-methylation.txt.gz \
        --regions r.bed

Note that the arguments are the same except for the --regions argument
giving the regions to test. This method can also work for e *X* pression.

Assumptions
===========

If you use the python module, these assumptions do not need to be met, you'll
just have to do some programming.

The command-line assumes that your methylation data looks like this:

```Shell
$ zless clustermodel/tests/example-methylation.txt.gz | cut -f 1-6 | head

probe   TF0 FF1 TM2 TF3 FM4
chr1:2041228    3.096   3.678   3.032   3.186   3.565
chr1:2041764    2.623   2.357   2.466   2.436   2.427
chr1:2043037    2.364   2.559   2.508   2.451   2.744
chr1:2043439    2.568   3.142   3.102   2.954   2.627
chr1:2043450    3.310   3.367   3.374   3.977   3.298
chr1:2043761    1.616   1.503   1.880   1.567   0.674
chr1:2043799    3.537   2.968   4.254   3.319   1.197
chr1:2043853    2.804   3.118   3.298   2.981   1.415
chr1:2045413    2.883   2.999   3.199   3.136   3.246
```
where rows are probes and columns are samples. The first column
should be in the form chrom:position and *must be sorted*
The columns of the methylation matrix must match the rows of the
covariates (here: `TF0`, `FF1`, `TM2`, `TF3`, `FM4`, ...):

```Shell
$ head -6 clustermodel/tests/example-covariates.txt

sample_id   disease gender  anumber
TF0 T   F   0
FF1 F   F   1
TM2 T   M   2
TF3 T   F   3
FM4 F   M   4
```

and the first column of the covariates file must match the first row
of the methylation data. (this should be pretty standard).
Any columns in the covariates file can be specified in the model.

So main points:

 1. data must be sorted
 2. methylation is n_probes * n_samples (+1 for index)
     a. first column is index of chrom:postion
     b. column headers are sample-ids matching first column from covariates
 3. covariates is n_samples * n_covariates (can have extra covarites even if
    not used in model.
     a. first column must match first row from methylation file


Methylation-eQTR
================
Methylation::expression Quantitative Trait Region. Generally, we perform a SNP-eQTL with
a model like:

    expression ~ genotype + age + gender ...

To find how gene expression is affected by genotype. We are always interested
in how methylation affects expression, but we often *lack power* to test each
*single methylation* probe against many expression probes. This software
enables testing the relation between expression and methylation in the context
of a mixed-effects model (or GEE, or liptak, or bumphunting).
The syntax looks like::

    python -m clustermodel \
           'methylation ~ disease + gender + (1|id) + (1|CpG)' \
           clustermodel/tests/example-covariates.txt \
           clustermodel/tests/example-methylation.txt.gz \
           --X clustermodel/tests/example-expression.txt.gz \
           --X-locs clustermodel/tests/example-expression-probe-locs.bed.gz \
           --X-dist 150000

where the first 3 arguments are as before: the model, the covariates, and the
methylation data.  The last 3 arguments are the expression info. The first
is the expression matrix with samples that match the covariates and the
methylation matrix. --X-locs is a BED file giving the location of each probe
and --X-dist gives the maximum distance between a methylation and expression
probe to be tested.  --X-locs and --X-dist are optional, but performing all
expression::methylation comparisons is likely going to take a long time,
despite the automatic parallelization.

The first few lines of output should look like::

    #chrom	start	end	coef	p	n_probes	model	method	Xstart	Xend	Xstrand	distance
    chr1	2043438	2043450	-0.05407	0.6498	2	methylation ~ A_24_P49214 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2116522	2116463	-	73072
    chr1	2043438	2043450	-0.01590	0.7547	2	methylation ~ A_33_P3410123 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	1920067	1920008	-	-123430
    chr1	2043438	2043450	-0.01726	0.8491	2	methylation ~ A_23_P51187 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2116703	2116762	+	-73253
    chr1	2043438	2043450	-0.12156	0.0962	2	methylation ~ A_33_P3359344 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2125138	2125079	-	81688
    chr1	2043438	2043450	0.030371	0.3701	2	methylation ~ A_33_P3410121 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	1915812	1915753	-	-127685
    chr1	2043438	2043450	-0.04481	0.2563	2	methylation ~ A_33_P3359354 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2116762	2116703	-	73312
    chr1	2043760	2043853	-0.05612	0.8589	3	methylation ~ A_24_P49214 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2116522	2116463	-	72669
    chr1	2043760	2043853	0.064542	0.6318	3	methylation ~ A_33_P3410123 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	1920067	1920008	-	-123752
    chr1	2043760	2043853	0.201022	0.4006	3	methylation ~ A_23_P51187 + disease + gender + (1 | id) + (1 | CpG)	mixed-model	2116703	2116762	+	-72850
    

Note that each probe is automatically inserted into the model. The distance
column is negative if the methylation region is upstream of the expression
probe and positive otherwise. 

See the example data in `clustermodel/tests/` for how to set up your own data.

Using Clustermodel With Bisulfite-Seq Data
==========================================

Most modes of use of `clustermodel` assume that the methylation data is a from a technology
that returns a proportion of methylated cells (450K, charm, etc), however, for sequencing
data, we know the proportion of methylation, but also the sequencing depth. To model this,
we use beta regression (similar to BiSeq from bioconductor). BiSeq and other
packages such as `bsseq` (also from bioconductor) use the sequencing depth as weights for
smoothing the methylation rate using loess and then model the smoothed methylation rates.
In `clustermodel`, we instead use the sequence depth to perform a weighted regression on
the original data such that **samples** with more supporting reads are weighted more heavily.
We then combine the p-values from each probe in a cluster (described below) using the z-score
method described in the BiSeq paper except that **sites** with more supporting reads are given
more weight.

To find clusters used above, we utilize the method described in http://www.ncbi.nlm.nih.gov/pubmed/23990415
and utilized throughout `clustermodel`. Again, this differs from `BiSeq` and `BSSeq` which
use peak-finding to delineate regions after asigning per-CpG p-values. The advantage of our approach is
that all clusters are reported so users can choose their own cutoff and multiple-testing correction is
more straight-forward.

An example invocation to use `clustermodel` would be:

```Shell

     python -m clustermodel \
            --png-path results/png/spaghetti \
            --counts --min-cluster-size 5 \
            --combine z-score --outlier-sds 3 \
             --betareg --weights counts.txt.gz \
            'methylation ~ case' \
             covariates.txt \
             methylation.txt.gz \
            > results/pvals.bed

Where `counts.txt.gz` has the same shape as `methylation.txt.gz` and the former has the read-depths
while the latter has the proportion of methylation for rows of sites and columns of samples.
The columns of those files will correspond to the rows in covariates.txt

The `counts` and `methylation` can be created by sending a list of Bismark or bwa-meth.py output files (tabulated methylation)
to `scripts/meth-matrix.py`.


```


INSTALLATION
============

You should have all of the R packages in scripts/check-R.R installed and should
be able to run that script without problems. You will have to do in R:

    install.packages(devtools)
    library(devtools)
    install_github("brentp/clustermodelr")

The rest can be installed with install.packages or biocLite().


Generating Correlated Data
==========================

see simulate.py

Coming Soon
===========

comparison
