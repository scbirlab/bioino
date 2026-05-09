# ── Outside GffFile (module level) ───────────────────────────────────────────

class _FeatureInterval:
    """Reference to a GffLine plus parameters for on-demand offset computation.

    Stores a reference to the original GffLine; no copy is made at build time.
    Call annotation_at(pos) to get a GffLine with offset and locus_tag set.
    """

    __slots__ = ('feature', 'locus_tag', 'offset_origin', 'sign')

    def __init__(self, feature: GffLine, locus_tag: str,
                 offset_origin: int, sign: float) -> None:
        self.feature = feature
        self.locus_tag = locus_tag
        self.offset_origin = offset_origin
        self.sign = sign  # +1.0 or -1.0

    def annotation_at(self, pos: int) -> GffLine:
        attrs = {
            **self.feature.attributes,
            'locus_tag': self.locus_tag,
            'offset': int((pos - self.offset_origin) * self.sign),
        }
        return replace(self.feature, attributes=attrs)


class _ChromosomeLookup:
    """Sorted-interval lookup for a single chromosome.

    Binary search gives O(log n) per query where n = number of intervals.
    Overlapping features are handled in O(log n + k), k = overlap depth.
    Uses a running-max array (_max_end) for early loop termination.

    Intervals must be added with non-decreasing start values.
    """

    __slots__ = ('_starts', '_ends', '_intervals', '_max_end')

    def __init__(self) -> None:
        self._starts: list = []
        self._ends: list = []
        self._intervals: list = []
        self._max_end: list = []  # _max_end[i] = max(_ends[0..i])

    def add(self, start: int, end: int, interval: '_FeatureInterval') -> None:
        self._starts.append(start)
        self._ends.append(end)
        self._intervals.append(interval)
        prev = self._max_end[-1] if self._max_end else 0
        self._max_end.append(max(prev, end))

    def at(self, pos: int) -> Tuple[GffLine, ...]:
        idx = bisect_right(self._starts, pos) - 1
        if idx < 0:
            return ()
        results = []
        i = idx
        while i >= 0 and self._max_end[i] >= pos:
            if self._ends[i] >= pos:
                results.append(self._intervals[i].annotation_at(pos))
            i -= 1
        return tuple(reversed(results))