"""Parsing columns of GFF3 files."""
from typing import Optional, Union
from dataclasses import asdict, dataclass, field

@dataclass
class GffColumns:
    """GFF-formatted columns.

    Attributes
    ==========
    seqid : str
        Name of chromosome.
    source : str
        Name of database or computer software source of annotation.
    feature : str
        Feature type, for example exon, gene, etc.
    start : str
        Start coordinate.
    end : int
        End coordinate.
    score : str, optional
        Score for feature. Default: ".".
    strand : str, optional
        Strandedness of feature. Either "+" or "-". Default: "+".
    phase : str or int, optional
        Location of first codon in feature relative to start. Default: ".".
        
    Methods
    =======
    __str__()
        Show the GFF-formatted columns.

    Examples
    ========
    >>> columns = "NC_000913.3   GenBank exon    1   100 .   +   .".split()
    >>> print(GffColumns(*columns))  # doctest: +NORMALIZE_WHITESPACE
    NC_000913.3 GenBank exon    1       100     .       +       .

    """
    seqid: str
    source: str
    feature: str
    start: Union[str, int]
    end: Union[str, int]
    score: Optional[Union[str, int]] = field(default=".")
    strand: Optional[str] = field(default="+")
    phase: Optional[Union[str, int]] = field(default=".")

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)

    def __str__(self) -> str:
        """Show the GFF-formatted columns."""
        return "\t".join(map(str, self.as_dict().values()))
    
    def as_dict(self) -> dict:
        return asdict(self)
