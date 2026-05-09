"""Genomic intervals utilities."""

from typing import TYPE_CHECKING, Any, Tuple
from bisect import bisect_right
from dataclasses import dataclass, field, replace

if TYPE_CHECKING:
    from .line import GffLine
else:
    GffLine = Any


@dataclass(frozen=True)
class FeatureInterval:
    """Reference to a GffLine.
    
    Includes parameters for on-demand offset computation.

    Stores a reference to the original GffLine. Call annotation_at(pos) 
    to get a GffLine with offset and locus_tag set.

    Examples
    ========
    >>> from bioino.gff.column import GffColumns
    >>> from bioino.gff.line import GffLine
    >>> from bioino.gff.intervals import FeatureInterval
    >>> line_plus = GffLine(GffColumns('chr1', 'src', 'gene', 10, 50, '.', '+', '.'), {'Name': 'gA'})
    >>> line_minus = GffLine(GffColumns('chr1', 'src', 'gene', 10, 50, '.', '-', '.'), {'Name': 'gB'})

    Offset is zero at the origin:

    >>> FeatureInterval(line_plus, 'gA', 10, 1.).annotation_at(10).attributes['offset']
    0

    Offset increases downstream on + strand:

    >>> FeatureInterval(line_plus, 'gA', 10, 1.).annotation_at(30).attributes['offset']
    20

    Offset is measured from end on - strand:

    >>> FeatureInterval(line_minus, 'gB', 50, -1.).annotation_at(50).attributes['offset']
    0
    >>> FeatureInterval(line_minus, 'gB', 50, -1.).annotation_at(30).attributes['offset']
    20

    locus_tag override is applied:

    >>> FeatureInterval(line_plus, '_down-gA', 10, 1.).annotation_at(60).attributes['locus_tag']
    '_down-gA'

    Original GffLine attributes are not mutated:

    >>> interval = FeatureInterval(line_plus, '_down-gA', 10, 1.)
    >>> _ = interval.annotation_at(60)
    >>> 'locus_tag' in line_plus.attributes
    False

    """
    feature: GffLine
    locus_tag: str
    offset_origin: int
    sign: float

    def annotation_at(
        self, 
        pos: int
    ) -> GffLine:
        attrs = self.feature.attributes | {
            "locus_tag": self.locus_tag,
            "offset": int((pos - self.offset_origin) * self.sign),
        }
        return replace(self.feature, attributes=attrs)


@dataclass
class ChromosomeLookup:
    """Sorted-interval lookup for a single chromosome.

    Binary search gives O(log n) per query where n = number of intervals.
    Overlapping features are handled in O(log n + k), k = overlap depth.
    Uses a running-max array (_max_end) for early loop termination.

    Intervals must be added with non-decreasing start values.

    Examples
    ========
    >>> from bioino.gff.column import GffColumns
    >>> from bioino.gff.line import GffLine
    >>> from bioino.gff.intervals import ChromosomeLookup, FeatureInterval
    >>> line = GffLine(GffColumns('chr1', 'src', 'gene', 10, 50, '.', '+', '.'), {'Name': 'gA'})
    >>> iv = FeatureInterval(line, 'gA', 10, 1.)
    >>> lookup = ChromosomeLookup()
    >>> lookup.add(10, 50, iv)

    Returns empty tuple outside all intervals:

    >>> lookup.at(5)
    ()
    >>> lookup.at(60)
    ()

    Returns annotation tuple inside interval:

    >>> result = lookup.at(30)
    >>> len(result)
    1
    >>> result[0].attributes['offset']
    20

    Overlapping intervals both returned:

    >>> line2 = GffLine(GffColumns('chr1', 'src', 'gene', 20, 70, '.', '+', '.'), {'Name': 'gB'})
    >>> lookup.add(20, 70, FeatureInterval(line2, 'gB', 20, 1.))
    >>> len(lookup.at(30))
    2
    
    """
    starts: list = field(default_factory=list)
    ends: list = field(default_factory=list)
    intervals: list = field(default_factory=list)
    max_end: list = field(default_factory=list)

    def add(
        self, 
        start: int, 
        end: int, 
        interval: FeatureInterval
    ) -> None:
        self.starts.append(start)
        self.ends.append(end)
        self.intervals.append(interval)
        prev = self.max_end[-1] if self.max_end else 0
        self.max_end.append(max(prev, end))

    def at(
        self, 
        pos: int
    ) -> Tuple[GffLine, ...]:
        idx = bisect_right(self.starts, pos) - 1
        if idx < 0:
            return ()
        results = []
        while idx >= 0 and self.max_end[idx] >= pos:
            if self.ends[idx] >= pos:
                results.append(self.intervals[idx].annotation_at(pos))
            idx -= 1
        return tuple(reversed(results))
