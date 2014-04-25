"""
make a weights and a methylation matrix from
bismark/bwa-meth style output
"""
import click
from toolshed import nopen
import heapq
import itertools as it
import os.path as op
import re

class Row(object):
    __slots__ = ('chrom', 'start', 'end', 'value', 'count', 'source')

    def __init__(self, toks, source=None):
        self.chrom = toks[0]
#               methylated=bismark[[5]],
#                            reads=bismark[[5]]+bismark[[6]])

        self.start, self.end = int(toks[1]), int(toks[2])
        self.count = int(toks[4]) + int(toks[5])
        self.value = "%.3f" % (float(toks[4]) / self.count)
        self.count = str(self.count)
        self.source = source

    def __cmp__(self, other):
        return cmp(self.chrom, other.chrom) or cmp(self.start, other.start)


def bed_merge(row_iterables, sources):
    assert len(sources) == len(row_iterables)

    for loc, cgs in it.groupby(heapq.merge(*row_iterables),
                            lambda cg: (cg.chrom, cg.start)):

        cgs = list(cgs)
        cg = cgs[0]
        present = dict((c.source, c) for c in cgs)

        # if a file doesn't have a record for here, just append 0
        values = [(present[s].value if s in present else 'NA') for s in sources]
        counts = [(present[s].count if s in present else '0') for s in sources]
        yield cg.chrom, cg.start, cg.end, values, counts

def gen_iterable(fname, source_from_fname):
    source = source_from_fname(fname)
    for i, toks in enumerate(x.rstrip("\r\n").split("\t") for x in nopen(fname)):
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        yield Row(toks, source)

@click.command()
@click.option('--prefix', help="output prefix", required=True)
@click.option("--name-re", default=r".+/(.+?).methylation.txt$",
            help="regexp to convert file name to sample name")
@click.option("--min-samples", default=2,
            help="skip sites where fewer files than this have coverage")
@click.argument('methylation-files', nargs=-1)
def main(prefix, name_re, min_samples, methylation_files):
    name_re = re.compile(r"%s" % name_re)
    if not prefix.endswith((".", "/")): prefix += "."
    fhm = nopen('{prefix}methylation.txt.gz'.format(prefix=prefix), 'w')
    fhc = nopen('{prefix}counts.txt.gz'.format(prefix=prefix), 'w')

    def source_from_fname(fname):
        try:
            return name_re.search(fname).groups(0)[0]
        except:
            return op.basename(fname)

    iterables = [gen_iterable(f, source_from_fname) for f in methylation_files]
    sources = [source_from_fname(f) for f in methylation_files]

    fmt = "{chrom}:{start}\t{vals}\n"
    fhm.write("probe\t%s" % "\t".join(sources) + "\n")
    fhc.write("probe\t%s" % "\t".join(sources) + "\n")
    for chrom, start, end, values, counts in bed_merge(iterables, sources):
        if sum(tryfloat(v) > 0 for v in values) < min_samples: continue
        vals = "\t".join(values)
        fhm.write(fmt.format(chrom=chrom, start=start, vals=vals))
        counts = "\t".join(counts)
        fhc.write(fmt.format(chrom=chrom, start=start, vals=counts))

def tryfloat(v):
    try: return float(v)
    except ValueError: return 0


if __name__ == "__main__":
    main()

