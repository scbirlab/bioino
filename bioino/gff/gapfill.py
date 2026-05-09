"""Utilities for making lookup tables from GFF files."""
from typing import TYPE_CHECKING, Any, Optional

from .intervals import ChromosomeLookup, FeatureInterval

if TYPE_CHECKING:
    from .line import GffLine
else:
    GffLine = Any

from . import DOWNSTREAM_PREFIX, NAME_ATTRIBUTE, UPSTREAM_PREFIX


def _gapfill_intervals(
    gff_line: GffLine,
    last_feature: Optional[GffLine] = None
) -> list:
    """Interval descriptors for the intergenic gap immediately before gff_line.

    The gap is split at its midpoint: the upstream half references last_feature
    (or gff_line itself when there is no preceding feature); the downstream
    half references gff_line. Returns a list of (start, end, FeatureInterval).
    """
    last_end = last_feature.columns.end if last_feature is not None else 0
    gap_start = last_end + 1
    gap_end = gff_line.columns.start - 1

    if gap_start > gap_end:
        return []

    gap_midpoint = gap_start + (gap_end - gap_start) // 2
    result = []

    # Upstream half
    anchor0 = last_feature if last_feature is not None else gff_line
    strand0 = anchor0.columns.strand
    if strand0 == "+":
        origin0, sign0 = anchor0.columns.start, 1.
        prefix0 = DOWNSTREAM_PREFIX if last_feature is not None else UPSTREAM_PREFIX
    else:
        origin0, sign0 = anchor0.columns.end, -1.
        prefix0 = UPSTREAM_PREFIX if last_feature is not None else DOWNSTREAM_PREFIX
    tag0 = prefix0 + anchor0.attributes.get(NAME_ATTRIBUTE, "")
    result.append((
        gap_start, 
        gap_midpoint,
        FeatureInterval(
            anchor0, 
            tag0, 
            origin0, 
            sign0,
        )))

    # Downstream half (may be empty if gap == 1)
    if (gap_midpoint + 1) <= gap_end:
        strand1 = gff_line.columns.strand
        if strand1 == "+":
            origin1, sign1, prefix1 = gff_line.columns.start, -1., UPSTREAM_PREFIX
        else:
            origin1, sign1, prefix1 = gff_line.columns.end, 1., DOWNSTREAM_PREFIX
        tag1 = prefix1 + gff_line.attributes.get(NAME_ATTRIBUTE, "")
        result.append((
            gap_midpoint + 1, 
            gap_end,
            FeatureInterval(
                gff_line, 
                tag1, 
                origin1, 
                sign1,
            ),
        ))
    return result
