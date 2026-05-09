"""Input and output functions and classes for GFF3 files.

Makes an attempt to conform to GFF3 but makes no guarantees.

"""

from typing import Dict, Iterable, Mapping, Optional, Tuple, Union

from collections import defaultdict
import csv
from dataclasses import asdict, dataclass, field, fields, replace
from io import TextIOWrapper
from itertools import chain
import sys

from carabiner import print_err
from carabiner.cast import cast
from tqdm.auto import tqdm

_GFF_COLNAMES = (
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
_GFF_FEATURE_BLOCKLIST = (
    "region", 
    "repeat_region",
)


def _cast_to_file_handle(file: Union[str, TextIOWrapper]) -> TextIOWrapper:
    if isinstance(file, TextIOWrapper):
        return file
    else:
        return cast(file, to=TextIOWrapper)


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
        return "\n".join(map(str, self.data))

    def write(
        self, 
        file: Optional[TextIOWrapper] = None
    ) -> None:
        """Write GFF-formatted line to file."""
        return print(str(self), file=file)


@dataclass
class GffColumns:
    """GFF-formatted columns.

    Attributes
    ----------
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
    -------
    __str__()
        Show the GFF-formatted columns.

    Examples
    --------
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
    

@dataclass
class GffFile:
    r"""Object for reading, writing, and manipulating GFF files.

    Attributes
    ----------
    lines : iterable of GffLines
        GFF data.
    metadata : GffMetadata, optional
        GFF header.
    lookup : bool, optional
        Whether to generate a lookup table based on the GFF file. 
        Default: False

    Methods
    -------
    as_dict()
        Produce an iterator of dictionaries.
    from_file()
        Read a GFF file.
    write()
        Write to a file.


    Examples
    --------
    >>> from io import StringIO
    >>> lines = ["##meta1 item1", 
    ...          "#meta2  item2  comment", 
    ...          "\t".join("test_seq    test_source gene    1   10  .   +   .   ID=test01;attr1=+".split()),
    ...          "\t".join("test_seq    test_source gene    9   100  .   +   .   Parent=test01;attr2=+".split())]
    >>> file = StringIO()
    >>> for line in lines:
    ...     print(line, file=file)
    >>> gff = GffFile.from_file(file)
    >>> gff.write()  # doctest: +NORMALIZE_WHITESPACE
    ##meta1 item1
    #meta2  item2  comment
    test_seq   test_source     gene    1       10      .       +       .       ID=test01;attr1=+
    test_seq   test_source     gene    9       100     .       +       .       Parent=test01;attr2=+
    
    """
    lines: Iterable[GffLine]
    metadata: Optional[Union[GffMetadata, Iterable[Union[Iterable, GffMetadatum]]]] = field(default_factory=list)
    lookup: Optional[bool] = field(default=False)
    _lookup: Dict[int, Tuple[GffLine]] = field(init=False, default_factory=dict)

    def __post_init__(self):
        if isinstance(self.metadata, Iterable):
            self.metadata = GffMetadata(self.metadata)
        if self.lookup:
            self.lines = tuple(self.lines)
            self._lookup = self._lookup_table()


    @staticmethod
    def _gapfill_table(
        gff_line: GffLine, 
        last_feature: Optional[GffLine] = None
    ) -> Dict[int, GffLine]:
        lookup_table = defaultdict(list)
        this_start = gff_line.columns.start
        if last_feature is None:
            last_end = 0
            intergenic0 = gff_line.copy()
        
        else:

            last_end = last_feature.columns.end
            intergenic0 = last_feature.copy()

        intergenic1 = gff_line.copy()

        gap_span = (this_start - 1) - (last_end + 1)
        gap_midpoint = last_end + 1 + gap_span // 2
        
        if intergenic0.columns.strand == "+":
            pre_mid_offset_start = intergenic0.columns.start
            pre_mid_sign = 1.
            pre_mid_prefix = "_down-" if last_feature is not None else "_up-" 
        else:
            pre_mid_offset_start = intergenic0.columns.end
            pre_mid_sign = -1.
            pre_mid_prefix = "_up-" if last_feature is not None else "_down-" 

        if intergenic1.columns.strand == "+":
            post_mid_offset_start = intergenic1.columns.start
            post_mid_sign = -1.
            post_mid_prefix = "_up-"
        else:
            post_mid_offset_start = intergenic1.columns.end
            post_mid_sign = 1.
            post_mid_prefix = "_down-"

        attr0 = intergenic0.attributes.copy() 
        attr0.update(dict(locus_tag=pre_mid_prefix + attr0["Name"]))
        attr1 = intergenic1.attributes.copy() 
        attr1.update(dict(locus_tag=post_mid_prefix + attr1["Name"]))
        # print(attr1)

        # fill in the gap
        for i in range(last_end + 1, gap_midpoint + 1):

            offset = (i - pre_mid_offset_start) * pre_mid_sign
            this_attr0 = attr0.copy()
            this_attr0.update(dict(offset=int(offset)))
            this_intergenic = replace(intergenic0, 
                                      attributes=this_attr0)

            lookup_table[i].append(this_intergenic)

        for i in range(gap_midpoint + 1, this_start):

            offset = (i - post_mid_offset_start) * post_mid_sign
            this_attr1 = attr1.copy()
            this_attr1.update(dict(offset=int(offset)))
            this_intergenic = replace(intergenic1, 
                                      attributes=this_attr1)

            lookup_table[i].append(this_intergenic)

        return lookup_table


    def _lookup_table(self) -> Dict[int, GffLine]:

        """Generate a lookup table for parent features in GFF.

        Results in a dictionary allowing lookup by chromosome location 
        to return feature annotations. Regions without annotation
        are automatically filled with references to upstream or 
        downstream features.

        Notes
        -----
        - Currently only works for single-chromosome files.
        - Only references parent features. Child features not yet indexed.
        - Will not work for GFFs with a single parent feature.
        - Ignores the following feature types: {}

        Returns
        -------
        dict
            Dictionary mapping chromosome location to feature annotation.

        """.format(", ".join(_GFF_FEATURE_BLOCKLIST))

        print_err("Building annotation lookup table.")

        lookup_table = defaultdict(list)

        last_feature = None 

        for gff_line in tqdm(self.lines):

            if (gff_line.columns.feature not in _GFF_FEATURE_BLOCKLIST and 
                "Name" in gff_line.attributes and
                "Parent" not in gff_line.attributes):
                    
                gap_table = self._gapfill_table(gff_line, 
                                                last_feature)
                
                lookup_table.update(gap_table)

                offset_start = (gff_line.columns.start 
                                if gff_line.columns.strand == "+" 
                                else gff_line.columns.end)

                for i in range(gff_line.columns.start, 
                               gff_line.columns.end + 1):

                    offset =  abs(i - offset_start)
                    this_gff_line = gff_line.copy()
                    this_gff_line.attributes["offset"] = offset

                    lookup_table[i].append(this_gff_line)

                last_feature = gff_line.copy()

        if last_feature.columns.strand == "+":
            last_offset_start = last_feature.columns.start
            last_sign = 1.
            last_prefix = "_down-"
        else:
            last_offset_start = last_feature.columns.end
            last_sign = -1.
            last_prefix = '_up-'

        attr = last_feature.attributes.copy()
        attr.update(dict(locus_tag=last_prefix + attr['Name']))

        for i in range(last_feature.columns.end, 
                       last_feature.columns.end + 1000):

            offset = (i - last_offset_start) * last_sign
            new_attributes = attr.copy() 
            new_attributes.update(dict(offset=int(offset)))
            last_feature = replace(last_feature, attributes=new_attributes)

            lookup_table[i].append(last_feature)

        missing_entries = [i for i in range(1, max(lookup_table) + 1) if len(lookup_table[i]) == 0]

        if len(missing_entries) > 0:

            raise AttributeError(f"Chromosome locations {', '.join(missing_entries)} are missing from lookup table.")

        return {key: tuple(value) for key, value in lookup_table.items()}


    def as_dict(self) -> Iterable[dict]:
        r"""Converts a `GffFile` to a stream of dictionaries.

        The resulting dicitonary from a `GffLine` has keys for the GFF columns 
        1-8 and each of the attributes from that line.

        Yields
        ------
        dict
            Dictionary corresponding to a GffLine.

        Examples
        --------
        >>> from io import StringIO
        >>> file = StringIO()
        >>> lines = ["TEST    test    gene    1       100     .       +       +  ID=test001;comment=Test".split(),
        ...          "TEST2    test2    gene    101       200     .       +       +  ID=test002;comment=Test2".split()]
        >>> for line in lines:
        ...     print('\t'.join(line), file=file)
        >>> list(GffFile.from_file(file).as_dict())  # doctest: +NORMALIZE_WHITESPACE
        [{'seqid': 'TEST', 'source': 'test', 'feature': 'gene', 'start': 1, 'end': 100, 'score': '.', 'strand': '+', 'phase': '+', 'ID': 'test001', 'comment': 'Test'}, {'seqid': 'TEST2', 'source': 'test2', 'feature': 'gene', 'start': 101, 'end': 200, 'score': '.', 'strand': '+', 'phase': '+', 'ID': 'test002', 'comment': 'Test2'}]
        
        """
        return (line.as_dict() for line in self.lines)

        
    def to_csv(
        self,
        file: TextIOWrapper = sys.stdout,
        write_metadata: bool = False,
        sep="","
    ) -> None:
        
        r"""Writes a `GffFile` to a delimited file.

        The resulting table has columns for the GFF columns 1-8 and 
        each of the unique attributes from the whole `GffFile`.
        
        Parameters
        ----------
        file : file-like
            File-like object which has a `write` method. Default: `sys.stdout`.
        write_metadata : bool, optional
            Whether to write metadata header. Default: False.
        sep : str, optional
            Delimiter between fields. Default: ','.

        Returns
        -------
        None

        Raises
        ------
        IOError
            If nothing to write.

        Examples
        --------
        >>> from io import StringIO
        >>> file = StringIO()
        >>> print('\t'.join('TEST    test    gene    1       100     .       +       +       ID=test001;comment=Test'.split()), 
        ...       file=file)
        >>> print('\t'.join('TEST    test    gene    121       120     .       +       -       ID=test001;tag=test_tag'.split()), 
        ...       file=file)
        >>> output = StringIO()
        >>> GffFile.from_file(file).to_csv(output)
        >>> output.seek(0)
        0
        >>> print("".join(output))  # doctest: +NORMALIZE_WHITESPACE
        seqid,source,feature,start,end,score,strand,phase,ID,comment,tag
        TEST,test,gene,1,100,.,+,+,test001,Test,
        TEST,test,gene,121,120,.,+,-,test001,,test_tag
        <BLANKLINE>
        >>> file = StringIO()
        >>> print('\t'.join('TEST    test    gene    1       100     .       +       +       ID=test001;comment=Test'.split()), 
        ...       file=file)
        >>> print('\t'.join('TEST    test    gene    121       120     .       +       -       ID=test001;tag=test_tag'.split()), 
        ...       file=file)
        >>> output = StringIO()
        >>> GffFile.from_file(file).to_csv(output, sep='\t')
        >>> output.seek(0)
        0
        >>> print("".join(output))  # doctest: +NORMALIZE_WHITESPACE
        seqid   source  feature start   end     score   strand  phase   ID      comment tag
        TEST    test    gene    1       100     .       +       +       test001 Test
        TEST    test    gene    121     120     .       +       -       test001         test_tag
        <BLANKLINE>

        """

        self.lines = tuple(self.lines)
        attribute_keys = set()
        main_cols = None

        print_err('Processing GFF attributes...')
        for i, gff_line in enumerate(tqdm(self.lines)):
            if i == 0:
                main_cols = [f.name for f in fields(gff_line.columns)]
            attribute_keys |= set(gff_line.attributes)

        if main_cols is None:
            raise IOError('GFF stream is empty.')

        if write_metadata:
            self.metadata.write(file=file)

        csv_fieldnames = list(chain(main_cols, sorted(attribute_keys)))
        writer = csv.DictWriter(
            file,
            fieldnames=csv_fieldnames,
            delimiter=sep,
        )
        writer.writeheader()
        for gff_line in self.as_dict():
            writer.writerow(gff_line)
        return None


    @staticmethod
    def _from_file(
        file: Union[str, TextIOWrapper]
    ) -> Iterable[Union[GffMetadata, GffLine]]:
        metadata = []
        metadata_shown = False
        with _cast_to_file_handle(file) as file:
            for line in file:
                line = line.strip()
                if line.startswith('#'):  # header
                    flag = 'constrained' if line.startswith('##') else 'free'
                    this_metadata = line.lstrip('#').lstrip().split('\t')
                    metadata.append(
                        GffMetadatum(
                            name=this_metadata[0], 
                            flag=flag, 
                            values=this_metadata[1:],
                        )
                    )
                elif len(line) > 0:  ## tab-delimited table
                    if not metadata_shown:
                        yield GffMetadata(metadata)
                    metadata_shown = True
                    data = line.split('\t')  ## Must be TAB otherwise columns 1-8 get messed up
                    try:
                        columns = GffColumns(*data[:8])
                    except TypeError:
                        print(data[:8], file=sys.stdout)
                        raise IOError(
                            '\n!!! ERROR: Probably corrupted file. '
                            'Here\'s the last line read:\n\n'
                            f'{line}\n\n'
                        )
                    attributes = '\t'.join(data[8:]) 
                    yield GffLine(columns, attributes)

    @classmethod
    def from_file(cls, 
                  file: Union[str, TextIOWrapper],
                  lookup: bool = False):
        """Stream records from a GFF file.

        Takes a file handle and instantiates a `GffFile` object.

        Parameters
        ----------
        file: TextIO
            File handle such as on generated by `open(f, mode='r')`.
        lookup: bool, optional
            Whether to create lookup table. Default: False.

        Returns
        -------
        GffFile

        Raises
        ------
        IOError
            If file is incorrectly formatted.

        """
        metadata = []
        file_parser = cls._from_file(file)
        for item in file_parser:
            if isinstance(item, GffMetadata):
                metadata = item
                break

        return cls(
            lines=(line for line in file_parser), 
            metadata=metadata,
            lookup=lookup,
        )
    
    def write(
        self, 
        file: Optional[TextIOWrapper] = None
    ) -> None:
        """Stream GFF records to a GFF file.
        
        Parameters
        ----------
        file : TextIO
            File handle such as on generated by `open(f, mode='w')`.

        """
        self.metadata.write(file=file)
        for line in self.lines:
            line.write(file=file)
        return None
