#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "click>=8.3.1",
#     "cyvcf2>=0.32.1",
# ]
# ///

from collections.abc import Iterator
from contextlib import contextmanager

import click
import cyvcf2


@contextmanager
def open_vcf(path) -> Iterator[cyvcf2.VCF]:
    """A context manager for opening a VCF file."""
    vcf = cyvcf2.VCF(path)
    try:
        yield vcf
    finally:
        vcf.close()


def variant_sites_are_equivalent(variant1, variant2):
    """Test if two variants are at the same site"""
    return variant1.CHROM == variant2.CHROM and variant1.POS == variant2.POS


def variant_alleles_are_equivalent(variant1, variant2):
    """Test if two variants represent equivalent alleles"""

    # by ID
    # return variant_sites_are_equivalent(variant1, variant2) and variant1.ID == variant2.ID

    # by REF/ALT
    return variant_sites_are_equivalent(variant1, variant2) and variant1.REF == variant2.REF and variant1.ALT == variant2.ALT


@click.command()
@click.argument("vcf_file", type=click.Path(exists=True))
@click.argument("variants_vcf_file", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    type=str,
    default="-",
    help="Where to write output file, defaults to standard out",
)
def cli(vcf_file, variants_vcf_file, output) -> None:
    """Normalise the variants in VCF_FILE to the defined list in VARIANTS_VCF_FILE.
    For each variant in VARIANTS_VCF_FILE, if VCF_FILE contains the variant
    then emit that record, otherwise emit an empty record at that site.
    """
    with open_vcf(vcf_file) as vcf, open_vcf(variants_vcf_file) as variants_vcf:
        writer = cyvcf2.Writer(output, vcf)
        writer.write_header()
        prev = None
        for variant in variants_vcf:
            if prev is not None and variant_sites_are_equivalent(variant, prev):
                v = prev
            else:
                v = next(vcf)
            if variant_alleles_are_equivalent(variant, v):
                writer.write_record(v)
                prev = None
            else:
                # add missing format fields
                # NOTE: assumes single sample
                variant_with_missing = writer.variant_from_string(f"{str(variant).rstrip()}\t.\t.")
                writer.write_record(variant_with_missing)
                prev = v


if __name__ == "__main__":
    cli()
