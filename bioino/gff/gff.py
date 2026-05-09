"""Input and output functions and classes for GFF3 files.

Makes an attempt to conform to GFF3 but makes no guarantees.

"""
from typing import Dict, Iterable, Optional, Tuple, Union

from collections import defaultdict
import csv
from dataclasses import dataclass, field, fields, replace
from io import TextIOWrapper
from itertools import chain
import sys

from carabiner import print_err
from tqdm.auto import tqdm

from . import DOWNSTREAM_PREFIX, NAME_ATTRIBUTE, UPSTREAM_PREFIX
from .column import GffColumns
from .gapfill import _gapfill_intervals
from .intervals import ChromosomeLookup, FeatureInterval
from .line import GffLine, GffMetadata, GffMetadatum
from .utils import _cast_to_file_handle

_GFF_FEATURE_BLOCKLIST: Tuple[str] = (
    "region", 
    "repeat_region",
)


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
    _lookup: Dict[str, ChromosomeLookup] = field(init=False, default_factory=dict)

    def __post_init__(self):
        if isinstance(self.metadata, Iterable):
            self.metadata = GffMetadata(self.metadata)
        if self.lookup:
            self.lines = tuple(self.lines)
            self._lookup = self._lookup_table()

    def _lookup_table(self) -> Dict[str, ChromosomeLookup]:
        """Build per-chromosome interval lookup.

        Results in a dictionary allowing lookup by chromosome location 
            to return feature annotations. Regions without annotation
            are automatically filled with references to upstream or 
            downstream features.

        Returns
        =======
        dict
            Dictionary mapping chromosome location to feature annotation.
            Intergenic regions are covered by half-gap
            intervals. Offsets and locus_tags are computed on demand at access.

        Notes
        =====
        - Now handles multi-chromosome GFFs.
        - Only references parent features. Child features not yet indexed.
        - Stores references to original GffLine objects; no per-position copies.
        - Ignores: {}
        """.format(', '.join(_GFF_FEATURE_BLOCKLIST))

        print_err("Building annotation lookup table.")
        tables: Dict[str, ChromosomeLookup] = {}
        last_feature: Dict[str, GffLine] = {}
        for gff_line in tqdm(self.lines):
            if any([
                gff_line.columns.feature in _GFF_FEATURE_BLOCKLIST,
                NAME_ATTRIBUTE not in gff_line.attributes,
                "Parent" in gff_line.attributes,
            ]):
                continue

            seqid = gff_line.columns.seqid
            if seqid not in tables:
                tables[seqid] = ChromosomeLookup()
            table = tables[seqid]
            prev = last_feature.get(seqid)

            for start, end, interval in _gapfill_intervals(gff_line, prev):
                table.add(start, end, interval)

            strand = gff_line.columns.strand
            origin = gff_line.columns.start if strand == '+' else gff_line.columns.end
            sign = 1. if strand == "+" else -1.
            locus_tag = gff_line.attributes.get(
                "locus_tag",
                gff_line.attributes.get(NAME_ATTRIBUTE, ''),
            )
            table.add(
                gff_line.columns.start,
                gff_line.columns.end,
                FeatureInterval(gff_line, locus_tag, origin, sign),
            )
            last_feature[seqid] = gff_line

        # Trailing 1000 bp past the last annotated feature on each chromosome
        for seqid, prev in last_feature.items():
            strand = prev.columns.strand
            if strand == '+':
                origin, sign, prefix = prev.columns.start, 1., DOWNSTREAM_PREFIX
            else:
                origin, sign, prefix = prev.columns.end, -1., UPSTREAM_PREFIX
            tag = prefix + prev.attributes.get('Name', '')
            tables[seqid].add(
                prev.columns.end + 1,
                prev.columns.end + 1000,
                FeatureInterval(prev, tag, origin, sign),
            )

        return tables

    def lookup_at(
        self, 
        seqid: str, 
        pos: int
    ) -> Tuple[GffLine, ...]:
        """Look up annotations at a chromosomal position.

        Parameters
        ==========
        seqid : str
            Chromosome identifier (i.e. GFF column 1).
        pos : int
            1-based chromosomal coordinate.

        Returns
        =======
        tuple of GffLine
            All features covering pos with offset and locus_tag computed.
            Empty tuple if pos falls outside all annotated intervals.

        # GffFile.lookup_at — goes in gff.py

        Examples
        ========
        >>> from io import StringIO
        >>> from bioino.gff.gff import GffFile
        >>> gff_text = '\\n'.join([
        ...     '\\t'.join(['chr1', 'src', 'gene', '10', '50', '.', '+', '.', 'ID=g1;Name=geneA']),
        ...     '\\t'.join(['chr1', 'src', 'gene', '100', '150', '.', '+', '.', 'ID=g2;Name=geneB']),
        ...     '\\t'.join(['chr2', 'src', 'gene', '20', '80', '.', '-', '.', 'ID=g3;Name=geneC']),
        ... ])
        >>> gff = GffFile.from_file(StringIO(gff_text), lookup=True)  # doctest: +ELLIPSIS
        ...

        Gene body, + strand — offset from start:

        >>> r = gff.lookup_at('chr1', 30)
        >>> r[0].attributes['locus_tag'], r[0].attributes['offset']
        ('geneA', 20)

        Intergenic, first half — attributed to upstream gene:

        >>> r = gff.lookup_at('chr1', 75)
        >>> r[0].attributes['locus_tag'], r[0].attributes['offset']
        ('_down-geneA', 65)

        Intergenic, second half — attributed to downstream gene:

        >>> r = gff.lookup_at('chr1', 76)
        >>> r[0].attributes['locus_tag'], r[0].attributes['offset']
        ('_up-geneB', 24)

        Trailing region past last gene:

        >>> r = gff.lookup_at('chr1', 200)
        >>> r[0].attributes['locus_tag'], r[0].attributes['offset']
        ('_down-geneB', 100)

        chr2 is independent — gene body on - strand, offset from end:

        >>> r = gff.lookup_at('chr2', 50)
        >>> r[0].attributes['locus_tag'], r[0].attributes['offset']
        ('geneC', 30)

        Unknown seqid returns empty tuple:

        >>> gff.lookup_at('chrX', 50)
        ()

        Position before any interval returns empty tuple:

        >>> gff.lookup_at('chr1', 0)
        ()
        
        """
        chrom = self._lookup.get(seqid)
        if chrom is None:
            return ()
        return chrom.at(pos)

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
        sep=","
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
        if len(self.metadata.data) > 0:
            self.metadata.write(file=file)
        for line in self.lines:
            line.write(file=file)
        return None
