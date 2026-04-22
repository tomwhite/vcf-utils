#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "click>=8.3.1",
#     "cyvcf2>=0.32.1",
#     "more-itertools>=10.8.0",
# ]
# ///

from collections.abc import Iterator
from contextlib import contextmanager

import click
import cyvcf2
from more_itertools import peekable


@contextmanager
def open_vcf(path) -> Iterator[cyvcf2.VCF]:
    """A context manager for opening a VCF file."""
    vcf = cyvcf2.VCF(path)
    try:
        yield vcf
    finally:
        vcf.close()


def variant_sites_equal(variant1, variant2):
    """Test if two variants have the same CHROM and POS"""
    return variant1.CHROM == variant2.CHROM and variant1.POS == variant2.POS

def sort_in_site(variants):
    it = peekable(variants)

    while True:
        if it.peek(None) is None:
            break
        # for variants at the same site, gather then merge
        v = next(it)
        vs = [v]
        try:
            while variant_sites_equal(v, it.peek()):
                vs.append(next(it))
        except StopIteration:
            pass

        if len(vs) == 1:
            yield vs[0]
        else:
            yield from sorted(vs, key=lambda x: (x.REF, tuple(x.ALT)))

@click.command()
@click.argument("vcf_file", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    type=str,
    default="-",
    help="Where to write output file, defaults to standard out",
)
def cli(vcf_file, output) -> None:
    """Sort variants at the same site by REF then ALT"""
    with open_vcf(vcf_file) as vcf:
        writer = cyvcf2.Writer(output, vcf)
        writer.write_header()

        for variant in sort_in_site(vcf):
            writer.write_record(variant)


if __name__ == "__main__":
    cli()