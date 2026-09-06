#!/usr/bin/env python
# (C) 2024, Tom Eulenfeld, MIT license
"""
anchorna-alimerge

Merge alignments with overlapping regions (anchors).

Anchors are cut from the start of all alignments except the first.
"""
from sugar import read


__version__ = '2024.09.05'


MAX = 100


def alimerge(alis, translate=True):
    """
    Merge alignments into single alignment

    :param alis: Input alignments
    :param translate: Wether to translate before identifying conserved
        regions, default is True.
    :returns: Merged alignment
    """
    if translate:
        talis = [ali.copy().translate(complete=True) for ali in alis]
    else:
        talis = [ali.copy() for ali in alis]
    ali1 = alis.pop(0)
    tali1 = talis.pop(0)
    nali = 0
    while len(alis) > 0:
        nali += 1
        ali2 = alis.pop(0)
        tali2 = talis.pop(0)
        ids = ali1.ids
        if set(ids) != set(ali2.ids):
            # TODO: handle this case
            raise ValueError('Not the same seqids in all files')
        n = None
        for i in range(1, min(min(len(seq) - seq.count('-')
                                  for seq in (tali1 + tali2)), MAX)):
            if all(tali1.d[id_].match(f'([^-]-*){{{i}}}$',
                                      gap=None).group() ==
                   tali2.d[id_].match(f'^(-*[^-]){{{i}}}',
                                      gap=None).group()
                   for id_ in ids):
                n = i
        if n is None:
            n = 0
        # TODO: log value of n
        matches = ali2.match(f'^(-*[^-]){{{3*n}}}')
        minoverlap = min(len(m.group()) for m in matches)
        for i in range(len(ali1)):
            m = matches[i]
            j = len(m.group())
            gaps = (j - minoverlap) * '-'
            ali1[i].data = ali1[i].data + gaps + ali2[i].data[j:]
        tali1 = tali2
    return ali1


def alimerge_io(fnames, output='-', fmt=None, output_fmt=None, translate=True):
    alis = [read(fname, fmt) for fname in fnames]
    seqs = alimerge(alis, translate=translate)
    if output_fmt is None and output == '-':
        output_fmt = 'fasta'
    seqs.write(output, fmt=output_fmt)


def test_simple():
    seqs1 = read()[:, :9]
    seqs2 = read()[:, 3:60]
    seqs3 = read()[:, 48:]
    alis = [seqs1, seqs2, seqs3]
    alimerge(alis)
    seqs = read()
    seqs[1].id = seqs[0].id
    ali1 = seqs[:1]
    ali1[0].data = ali1[0].data + seqs[1].data
    ali2 = alimerge([seqs[:1], seqs[1:]])
    assert ali1 == ali2


def test_hardcore():
    """Run tests with anchorna-alicat --test"""
    from sugar.tests.util import _changetmpdir
    from anchorna.cli import run_cmdline as run
    with _changetmpdir():
        print('Run AnchoRNA tutorial subset')
        run('create --tutorial-subset'.split())
        run('go a.gff'.split())
        run('print a.gff'.split())
        print('\nCut out stuff and check')
        run('cutout a.gff ATG A0 -o ali1.fasta'.split())
        run('cutout a.gff A0 A1 -o ali2.fasta'.split())
        run('cutout a.gff A1 A2 -o ali3.fasta'.split())
        alimerge_io(('ali1.fasta', 'ali2.fasta', 'ali3.fasta'), output='out.fasta')
        run('cutout a.gff ATG A2 -o alitest.fasta'.split())
        out1 = read('out.fasta')
        out2 = read('alitest.fasta')
        assert out1 == out2
        print(out1)
        print('OK')


def main(args=None):
    import argparse
    import sys

    if '--test' in sys.argv:
        test_simple()
        test_hardcore()
        sys.exit()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    p.add_argument('fnames', nargs='+', help='alignment files (or seqs files), sorted from left to right')
    p.add_argument('-o', '--output', help='output file name (default: pipe)', default='-')
    p.add_argument('--fmt', help='alignment format (default: auto-detect)')
    p.add_argument('--output-fmt', help='output format (default: auto-detect)')
    p.add_argument('--no-translate', action='store_false', dest='translate', help='turn of translation')
    p.add_argument('--test', default=argparse.SUPPRESS, help='run tests, all other arguments are ignored')
    args = p.parse_args(args)
    alimerge_io(**vars(args))
