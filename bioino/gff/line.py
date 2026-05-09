"""Parsing lines of GFF3 files.

"""
from typing import Dict, Iterable, Mapping, Optional, Tuple, Union

from dataclasses import dataclass, field, fields, replace
from io import TextIOWrapper
import sys

from .column import GffColumns

_GFF_COLNAMES: Tuple[str] = (
    "seqid", 
    "source", 
    "feature", 
    "start", 
    "end", 
    "score", 
    "strand", 
    "phase", 
    "attribute",
)

@dataclass
class GffMetadatum:
    """GFF-formatted metadata line.

    Attributes
    ----------
    name : str
        Name of metadatum.
    flag : str, optional
        "constrained" or "free", depending on whether it conforms to GFF. Default: "free".
    values : tuple, optional
        Tuple of values corresponding to `name`. Default: zero-length tuple.
        
    Methods
    -------
    __str__()
        Show the GFF-formatted metadata.
    write()
        Write GFF-formatted line to file.
        

    Examples
    --------
    >>> print(GffMetadatum("Meta_name", "free", ("meta_value1", "meta_value2")))  # doctest: +NORMALIZE_WHITESPACE
    #Meta_name  meta_value1     meta_value2
    >>> print(GffMetadatum("Meta_name", "constrained", ("meta_value1", "meta_value2")))  # doctest: +NORMALIZE_WHITESPACE
    ##Meta_name meta_value1     meta_value2

    """
    name: str
    flag: str = field(default="constrained")
    values: Tuple = field(default_factory=tuple)

    def __post_init__(self):
        if self.flag not in ["free", "constrained"]:
            raise ValueError("GffMetadatum.flag must be one of ['free', 'constrained']].")
        
    def __str__(self) -> str:
        """Show the GFF-formatted metadata."""
        prefix = "##" if self.flag == "constrained" else "#"
        suffix = "\t".join(map(str, self.values))
        return f"{prefix}{self.name}\t{suffix}"
    
    def write(
        self, 
        file: TextIOWrapper = sys.stdout
    ) -> None:
        """Write GFF-formatted line to file."""
        return print(str(self), file=file)
        

@dataclass
class GffMetadata:
    """GFF-formatted metadata.

    Attributes
    ----------
    metadata : Tuple[GffMetadatum]
        Tuple of metadata lines.

    Methods
    -------
    __str__()
        Show the GFF-formatted metadata.
    write()
        Write GFF-formatted line to file.

    Examples
    --------
    >>> metadata = [("meta1", "constrained", ("item1", )), 
    ...             ("meta2", "free", ("item2", "comment"))]
    >>> metadata = GffMetadata(metadata)
    >>> print(metadata)  # doctest: +NORMALIZE_WHITESPACE
    ##meta1 item1
    #meta2 item2    comment

    """
    data: Iterable[Union[GffMetadatum, Iterable]]

    def __post_init__(self):
        new_metadata = []
        for item in self.data:
            if isinstance(item, GffMetadatum):
                new_metadatum = item
            elif isinstance(item, Iterable) and not isinstance(item, str):
                new_metadatum = GffMetadatum(*item)
            else:
                raise ValueError(f"{item} of type {type(item)} cannot be converted to GffMetadatum.")
            new_metadata.append(new_metadatum)
        self.data = tuple(new_metadata)
    
    def __str__(self) -> str:
        """Show the GFF-formatted metadata."""
        if len(self.data) > 0:
            return "\n".join(map(str, self.data))
        else:
            return ""

    def write(
        self, 
        file: Optional[TextIOWrapper] = None
    ) -> None:
        """Write GFF-formatted line to file."""
        return print(str(self), file=file)


@dataclass
class GffLine:
    """Named tuple which gives a GFF-formatted line when printed.

    Attributes
    ----------
    metadata : tuple
        Tuple of GffMetadata from the original file.
    columns : GffColumns
        Representation of columns 1-8.
    attributes : dict
        Dictionary mapping attribute keys to values.

    Methods
    -------
    copy()
        Make a copy.
    __str__()
        Show the GFF-formatted line.

    Examples
    --------
    >>> metadata = [("meta1", "constrained", {"item1": []}), 
    ...             ("meta2", "free", {"item2": ["comment"]})]
    >>> columns = ["test_seq", "test_source", "gene", 1, 10]
    >>> gff_line = GffLine(columns, 
    ...                    attributes={"ID": "test01", "attr1": "+"})  
    >>> print(gff_line)  # doctest: +NORMALIZE_WHITESPACE
    test_seq        test_source     gene    1       10      .       +       .       ID=test01;attr1=+

    """

    columns : Union[GffColumns, Iterable]
    attributes : Optional[dict] = field(default_factory=dict)

    @staticmethod
    def _get_gff_attributes(x: str) -> Dict[str, str]:
        splits_on_equal_sign = [item.split(";") for item in x.split("=")]
        attributes = (item[-1] for item in splits_on_equal_sign)
        values = (item[0] for item in splits_on_equal_sign[1:])
        return dict(zip(attributes, values))
    
    def __post_init__(self):
        if isinstance(self.columns, Iterable):
            self.columns = GffColumns(*self.columns)
        if isinstance(self.attributes, str):
            self.attributes = self._get_gff_attributes(self.attributes)
    
    def __str__(self) -> str:
        """Show the GFF-formatted line."""
        _attributes = ";".join(f"{key}={val}" for key, val in self.attributes.items())
        return str(self.columns) + "\t" + _attributes
    
    def as_dict(self) -> dict:
        """Convert to dictionary."""
        d = self.columns.as_dict()
        d.update(self.attributes)
        return d

    def copy(self):
        """Make a copy."""
        return replace(self)
    
    def write(
        self,
        file: Optional[TextIOWrapper] = None
    ) -> None:
        """Write GFF-formatted line to file."""
        return print(str(self), file=file)
    

    @classmethod
    def from_dict(
        cls,
        d: Mapping
    ):
        """Converts a dictionary object to a GFFLine.

        The input dictionary must at least have keys corresponding to the GFF
        columns 1-8 and optionally additionl keys to put into the attributes
        column.
        
        Parameters
        ----------
        d : dict
            Dictionary to convert.

        Returns
        -------
        GffLine
            Object representing a GFF line.

        Examples
        --------
        >>> d = dict(seqid="TEST", source="test", 
        ...          feature="gene", start=1, 
        ...          end=100, score=".", 
        ...          strand="+", phase="+")
        >>> print(GffLine.from_dict(d)) # doctest: +NORMALIZE_WHITESPACE
        TEST        test    gene    1       100     .       +       +
        >>> d.update(dict(ID="test001", comment="This is a test"))
        >>> GffLine.from_dict(d).write() # doctest: +NORMALIZE_WHITESPACE
        TEST    test    gene    1       100     .       +       +       ID=test001;comment=This is a test

        """
        _fields = _GFF_COLNAMES #[f.name for f in fields(cls)]
        columns = GffColumns(**{
            key: value for key, value in d.items() 
            if key in _fields
        })
        attributes = {
            key: d[key] for key in sorted(d) 
            if key not in _fields
        }
        return cls(columns, attributes)