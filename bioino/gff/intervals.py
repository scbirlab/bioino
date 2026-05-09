"""Genomic intervals utilities."""
from typing import TYPE_CHECKING, Any, 
import bisect
from dataclasses import dataclass

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
            return None
        results = []
        while idx >= 0 and self.max_end[idx] >= pos:
            if self._ends[idx] >= pos:
                results.append(self.intervals[idx].annotation_at(pos))
            idx -= 1
        return tuple(reversed(results))
