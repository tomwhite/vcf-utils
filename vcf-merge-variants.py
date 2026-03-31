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

from typing import Any, Callable

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


def cmp_variants(variant1, variant2):
    if variant1.CHROM != variant2.CHROM:
        raise ValueError("Cannot compare variants with different CHROM")
    return int(variant1.POS > variant2.POS) - int(variant1.POS < variant2.POS) 


def variant_is_not_after(variant1, variant2):
    """Test if variant 1 is not after variant 2 along the genome"""
    return variant1.CHROM == variant2.CHROM and variant1.POS <= variant2.POS


def variant_alleles_are_equivalent(variant1, variant2):
    """Test if two variants represent equivalent alleles"""

    # by ID
    # return variant1.CHROM == variant2.CHROM and variant1.POS == variant2.POS and variant1.ID == variant2.ID

    # by REF/ALT
    return variant1.CHROM == variant2.CHROM and variant1.POS == variant2.POS and variant1.REF == variant2.REF and variant1.ALT == variant2.ALT


def merge_lists(l1: list, l2: list, *, key: Callable[[Any], Any] = lambda x: x) -> list:
    """Merge two lists preserving the relative order of elements from both inputs.

    Returns a list containing all elements from l1 and l2, deduplicated by key,
    ordered consistently with both inputs. Original items (not key values) are returned;
    when the same key appears in both lists, the item from l1 is kept.

    Args:
        l1: First ordered list.
        l2: Second ordered list.
        key: Function extracting a hashable identity from each element. Defaults to
             the element itself. Two elements with the same key are considered equal.

    Raises ValueError if the ordering constraints from l1 and l2 conflict,
    or if either input list contains duplicate elements (by key).
    """
    for lst, name in ((l1, "l1"), (l2, "l2")):
        keys = [key(item) for item in lst]
        if len(keys) != len(set(keys)):
            raise ValueError(f"Input {name} contains duplicate elements")

    # Assign a stable rank to each element keyed by key(item).
    # When the same key appears in both lists, the item from l1 is kept.
    rank: dict[Any, int] = {}
    items: dict[Any, Any] = {}  # key -> original item
    for item in l1 + l2:
        k = key(item)
        if k not in rank:
            rank[k] = len(rank)
            items[k] = item

    all_keys = list(rank.keys())
    graph: dict[Any, set] = {k: set() for k in all_keys}
    in_degree: dict[Any, int] = {k: 0 for k in all_keys}

    for lst in (l1, l2):
        for a, b in zip(lst, lst[1:]):
            ka, kb = key(a), key(b)
            if kb not in graph[ka]:
                graph[ka].add(kb)
                in_degree[kb] += 1

    queue = sorted([k for k in all_keys if in_degree[k] == 0], key=rank.__getitem__)
    result_keys: list = []

    while queue:
        k = queue.pop(0)
        result_keys.append(k)
        newly_free = []
        for successor in graph[k]:
            in_degree[successor] -= 1
            if in_degree[successor] == 0:
                newly_free.append(successor)
        queue = sorted(queue + newly_free, key=rank.__getitem__)

    if len(result_keys) < len(all_keys):
        remaining = [items[k] for k in all_keys if k not in set(result_keys)]
        raise ValueError(
            f"Cannot merge lists: ordering conflict detected among elements {remaining}"
        )
    return [items[k] for k in result_keys]


def merge(variants1, variants2):
    it1, it2 = peekable(variants1), peekable(variants2)

    while True:
        v1 = it1.peek(None)
        if v1 is None:
            yield from it2
            break
        v2 = it2.peek(None)
        if v2 is None:
            yield from it1
            break
        c = cmp_variants(v1, v2)
        if c < 0:
            yield next(it1)
        elif c > 0:
            yield next(it2)
        else:
            # for variants at the same site, gather then merge
            v1 = next(it1)
            v1s = [v1]
            try:
                while cmp_variants(v1, it1.peek()) == 0:
                    v1s.append(next(it1))
            except StopIteration:
                pass

            v2 = next(it2)
            v2s = [v2]
            try:
                while cmp_variants(v2, it2.peek()) == 0:
                    v2s.append(next(it2))
            except StopIteration:
                pass

            merged = merge_lists(v1s, v2s, key=lambda x: (x.REF, tuple(x.ALT)))
            yield from merged


@click.command()
@click.argument("vcf_file1", type=click.Path(exists=True))
@click.argument("vcf_file2", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    type=str,
    default="-",
    help="Where to write output file, defaults to standard out",
)
def cli(vcf_file1, vcf_file2, output) -> None:
    """Merge the variants in VCF_FILE1 with VCF_FILE2"""
    with open_vcf(vcf_file1) as vcf1, open_vcf(vcf_file2) as vcf2:
        writer = cyvcf2.Writer(output, vcf1)
        writer.write_header()

        for variant in merge(vcf1, vcf2):
            writer.write_record(variant)


if __name__ == "__main__":
    cli()
