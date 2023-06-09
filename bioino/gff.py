"""Input and output functions and classes for GFF3 files.

Makes an attempt to conform to GFF3 but makes no guarantees.

"""

from typing import TextIO, Union
from collections import defaultdict, namedtuple
from collections.abc import Generator, Iterable, Mapping, Sequence
import csv
import io
import sys

_GFF_COLNAMES = ('seqid', 'source', 'feature', 
                'start', 'end', 'score', 
                'strand', 'phase', 'attribute')
_GFF_FEATURE_BLOCKLIST = ('region', 'repeat_region')

_GffLine = namedtuple('_GffLine', 'metadata columns attributes')
_GffMetadatum = namedtuple('GffMetadatum', 'name flag values')
_GffColumns = namedtuple('GffColumns', _GFF_COLNAMES[:-1])

class GffLine(_GffLine):

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
    copy
        Make a copy.
    __str__
        Show the GFF-formatted line.

    Examples
    --------
    >>> metadata = [("meta1", "constrained", {"item1": []}), 
    ...             ("meta2", "free", {"item2": ["comment"]})]
    >>> seqid, source_id = "test_seq", "test_source"
    >>> print(GffLine(metadata, 
    ...               GffColumns(seqid, source_id, "gene", 1, 10, ".", "+", "."), 
    ...               {"ID": "test01", "attr1": "+"}))  # doctest: +NORMALIZE_WHITESPACE
    test_seq        test_source     gene    1       10      .       +       .       ID=test01;attr1=+
    <BLANKLINE>

    """
    
    def copy(self) -> Sequence:

        """Make a copy."""

        return super()._make([self.metadata, self.columns, self.attributes.copy()])
    
    def __str__(self) -> str:

        """Show the GFF-formatted line."""

        _attributes = ';'.join([f'{key}={val}' for key, val in self.attributes.items()])

        return (str(self.columns) + '\t' + _attributes + '\n')


class GffMetadatum(_GffMetadatum):

    """Named tuple which gives a GFF-formatted metadata 
    line when printed.

    Attributes
    ----------
    name : str
        Name of metadatum.
    flag : str
        'constrained' or 'free', depending on whether it conformas to GFF.
    values : tuple
        Tuple of values corresponding to `name`.
        
    Methods
    -------
    __str__()
        Show the GFF-formatted metadata.

    Examples
    --------
    >>> print(GffMetadatum('Meta_name', 'free', ('meta_value1', 'meta_value2')))  # doctest: +NORMALIZE_WHITESPACE
    #Meta_name  meta_value1     meta_value2
    >>> print(GffMetadatum('Meta_name', 'constrained', ('meta_value1', 'meta_value2')))  # doctest: +NORMALIZE_WHITESPACE
    ##Meta_name meta_value1     meta_value2

    """

    def __str__(self) -> str:

        """Show the GFF-formatted metadata."""

        prefix = ('##' if self.flag == 'constrained' 
                  else '#')

        return (prefix + self.name + '\t' + 
                '\t'.join(map(str, self.values)))
        

class GffColumns(_GffColumns):

    """Named tuple which gives a GFF-formatted metadata 
    line when printed.

    See GFF-format documentation for full description of
    attributes.

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
    score : str
        Score for feature.
    strand : str
        Strandedness of feature. Either '+' or '-'.
    phase : str or int
        Location of first codon in feature relative to start.
        
    Methods
    -------
    __str__()
        Show the GFF-formatted columns.

    Examples
    --------
    >>> print(GffColumns(*"NC_000913.3   GenBank exon    1   100 .   +   .".split()))  # doctest: +NORMALIZE_WHITESPACE
    NC_000913.3 GenBank exon    1       100     .       +       .

    """
    
    def __str__(self) -> str:

        """Show the GFF-formatted columns."""

        return '\t'.join(map(str, self))


def _get_gff_attributes(x: str) -> Mapping:

    splits_on_equal_sign = tuple(item.split(';') for item in x.split('='))
    
    attributes = (item[-1] for item in splits_on_equal_sign)
    values = (item[0] for item in splits_on_equal_sign[1:])
    
    return dict(zip(attributes, values))


def read_gff(file: Union[str, TextIO]) -> Generator[GffLine]:

    """Stream records from a GFF file.

    Takes a file handle and yields tuple of (metadata, columns, attributes) 
    for each record.

    Parameters
    ----------
    file: TextIO
        File handle such as on generated by `open(f, mode='r')`.

    Yields
    ------
    GffLine
        Named tuple containing metadata, columns, attributes fields.
    metadata : list
        GffMetadatum objects in the same order as the GFF header.
    columns : GffColumns
        Named tuple of columns 1 to 8 of the GFF.
    attributes : dict
        Dictionary of attributes found in column 9 of the GFF.

    Examples
    --------
    >>> list(read_gff('''
    ...     ##meta1 item1
    ...     #meta2  item2   comment
    ...     test_seq    test_source gene    1   10  .   +   .   ID=test01;attr1=+
    ...     test_seq    test_source gene    9   100  .   +   .   Parent=test01;attr2=+
    ...     '''))
    [GffLine(metadata=[GffMetadatum(name='meta1 item1', flag='constrained', values=()), GffMetadatum(name='meta2  item2   comment', flag='free', values=())], columns=GffColumns(seqid='test_seq', source='test_source', feature='gene', start=1, end=10, score='.', strand='+', phase='.'), attributes={'ID': 'test01', 'attr1': '+'}), GffLine(metadata=[GffMetadatum(name='meta1 item1', flag='constrained', values=()), GffMetadatum(name='meta2  item2   comment', flag='free', values=())], columns=GffColumns(seqid='test_seq', source='test_source', feature='gene', start=9, end=100, score='.', strand='+', phase='.'), attributes={'Parent': 'test01', 'attr2': '+'})]

    
    """

    if isinstance(file, str):
        file = io.StringIO(file)

    metadata = []

    for line in file:

        line = line.strip()

        if line.startswith('#'):  # header

            flag = 'constrained' if line.startswith('##') else 'free'
            this_metadata = line.lstrip('#').lstrip().split('\t')
            
            metadata.append(GffMetadatum(name=this_metadata[0], 
                                         flag=flag, 
                                         values=tuple(this_metadata[1:])))

        elif len(line) > 0:  # tab-delimited table

            data = line.split()
            
            try:
                
                columns = GffColumns(*data[:8])
            
            except TypeError as e:
                
                sys.stderr.write('\n!!! ERROR: Probably corrupted file. '
                                 'Here\'s the last line read:'
                                 '\n\n{}\n\n'.format(line))
                
                raise e
            
            columns = columns._replace(start=int(columns.start),
                                       end=int(columns.end))
            col9 = '\t'.join(data[8:])

            attributes = _get_gff_attributes(col9)

            yield GffLine(metadata, columns, attributes)


def write_gff(stream: Sequence[GffLine],
              file: TextIO = sys.stdout,
              metadata: bool = False) -> None:
    
    """Stream GFF records to a GFF file.

    Takes an iterable of FastaSequence and writes them to the given file.

    Parameters
    ----------
    stream : Sequence
        Iterable of GffLine objects.
    file : TextIO
        File handle such as on generated by `open(f, mode='w')`.
    metadata : bool
        Whether to write a metadata header. 

    Examples
    --------
    >>> metadata = [GffMetadatum("meta1", "constrained", {"item1": []}), 
    ...             GffMetadatum("meta2", "free", {"item2": ["comment"]})]
    >>> seqid, source_id = "test_seq", "test_source"
    >>> gff_stream = [GffLine(metadata, 
    ...                       GffColumns(seqid, source_id, "gene", 1, 10, ".", "+", "."), 
    ...                       {"ID": "test01", "attr1": "+"}),
    ...               GffLine(metadata, 
    ...                       GffColumns(seqid, source_id, "CDS", 9, 100, ".", "+", "."), 
    ...                       {"Parent": "test01", "attr2": "-"})]
    >>> write_gff(gff_stream)  # doctest: +NORMALIZE_WHITESPACE +SKIP
    test_seq        test_source     gene    1       10      .       +       .       ID=test01;attr1=+
    test_seq        test_source     CDS     9       100     .       +       .       Parent=test01;attr2=-
    >>> write_gff(gff_stream, metadata=True)  # doctest: +NORMALIZE_WHITESPACE +SKIP
    ##meta1 item1
    #meta2  item2
    test_seq        test_source     gene    1       10      .       +       .       ID=test01;attr1=+
    test_seq        test_source     CDS     9       100     .       +       .       Parent=test01;attr2=-

    """

    for i, line in enumerate(stream):

        if metadata and i == 0:

            file.write('\n'.join(map(str, line.metadata)) + '\n')

        file.write(str(line))

    return None


def _gapfill_table(gff_line: GffLine, 
                   last_feature: GffLine = None) -> Mapping:

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
    
    if intergenic0.columns.strand == '+':
        pre_mid_offset_start = intergenic0.columns.start
        pre_mid_sign = 1.
        pre_mid_prefix = '_down-' if last_feature is not None else '_up-' 
    else:
        pre_mid_offset_start = intergenic0.columns.end
        pre_mid_sign = -1.
        pre_mid_prefix = '_up-' if last_feature is not None else '_down-' 

    if intergenic1.columns.strand == '+':
        post_mid_offset_start = intergenic1.columns.start
        post_mid_sign = -1.
        post_mid_prefix = '_up-'
    else:
        post_mid_offset_start = intergenic1.columns.end
        post_mid_sign = 1.
        post_mid_prefix = '_down-'

    attr0 = intergenic0.attributes.copy() 
    attr0 |= dict(locus_tag=pre_mid_prefix + attr0['Name'])
    attr1 = intergenic1.attributes.copy() 
    attr1 |= dict(locus_tag=post_mid_prefix + attr1['Name'])
    # print(attr1)

    # fill in the gap
    for i in range(last_end + 1, gap_midpoint + 1):

        offset = (i - pre_mid_offset_start) * pre_mid_sign
        # attr0.update(dict(offset=int(offset)))
        this_intergenic = intergenic0._replace(attributes=(attr0 | dict(offset=int(offset))))

        lookup_table[i].append(this_intergenic)

    for i in range(gap_midpoint + 1, this_start):

        offset = (i - post_mid_offset_start) * post_mid_sign
        # attr1.update(dict(offset=int(offset)))
        this_intergenic = intergenic1._replace(attributes=(attr1 | dict(offset=int(offset))))

        lookup_table[i].append(this_intergenic)

    return lookup_table


def lookup_table(gff_stream: Iterable[GffLine]) -> Mapping:

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

    Parameters
    ----------
    gff_stream : Iterable[GffLine]
        An iterable of GffLine objects.

    Returns
    -------
    dict
        Dictionary mapping chromosome location to feature annotation.

    """.format(', '.join(_GFF_FEATURE_BLOCKLIST))

    lookup_table = defaultdict(list)

    last_feature = None 

    for gff_line in gff_stream:

        if (gff_line.columns.feature not in _GFF_FEATURE_BLOCKLIST and 
            'Name' in gff_line.attributes and
            'Parent' not in gff_line.attributes):
                
            gap_table = _gapfill_table(gff_line, 
                                       last_feature)
            
            lookup_table |= gap_table

            offset_start = (gff_line.columns.start 
                            if gff_line.columns.strand == '+' 
                            else gff_line.columns.end)

            for i in range(gff_line.columns.start, 
                           gff_line.columns.end + 1):

                offset =  abs(i - offset_start)
                this_gff_line = gff_line.copy()
                this_gff_line.attributes['offset'] = offset

                lookup_table[i].append(this_gff_line)

            last_feature = gff_line.copy()

    if last_feature.columns.strand == '+':
        last_offset_start = last_feature.columns.start
        last_sign = 1.
        last_prefix = '_down-'
    else:
        last_offset_start = last_feature.columns.end
        last_sign = -1.
        last_prefix = '_up-'

    attr = last_feature.attributes.copy() 
    attr |= dict(locus_tag=last_prefix + attr['Name'])

    for i in range(last_feature.columns.end, 
                   last_feature.columns.end + 1000):

        offset = (i - last_offset_start) * last_sign
        last_feature = last_feature._replace(attributes=(attr.copy() | dict(offset=int(offset))))

        lookup_table[i].append(last_feature)

    assert all(len(lookup_table[i]) > 0 for i in range(1, max(lookup_table) + 1))

    return lookup_table


def dict2gff(x: Mapping, 
             metadata: GffMetadatum = GffMetadatum('','','')) -> GffLine:
    
    """Converts a dictionary object to a GFFLine.

    The input dictionary must at least have keys corresponding to the GFF
    columns 1-8 and optionally additionl keys to put into the attributes
    column.
    
    Parameters
    ----------
    x : dict
        Dictionary to convert.
    metadata : GffMetadatum, optional
        Metadata lines to add to GFF.

    Returns
    -------
    GffLine
        Object representing a GFF line.

    Examples
    --------
    >>> d = dict(seqid='TEST', source='test', feature='gene', start=1, end=100, score='.', strand='+', phase='+')
    >>> print(dict2gff(d))  # doctest: +NORMALIZE_WHITESPACE
    TEST        test    gene    1       100     .       +       +
    <BLANKLINE>
    >>> print(dict2gff(d | dict(ID='test001', comment='This is a test')))  # doctest: +NORMALIZE_WHITESPACE
    TEST    test    gene    1       100     .       +       +       ID=test001;comment=This is a test
    <BLANKLINE>

    """

    columns = GffColumns(**{key: value for key, value in x.items() 
                            if key in GffColumns._fields})
    attributes = {key: x[key] for key in sorted(x) 
                  if key not in GffColumns._fields}

    return GffLine(metadata, columns, attributes)


def gff2dict(gff_stream: Iterable[GffLine]) -> Generator[Mapping]:

    """Converts a stream of GFFLines to a stream of dictionaries.

    The resulting dicitonary from a GffLine has keys for the GFF columns 
    1-8 and each of the attributes from that line.
    
    Parameters
    ----------
    x : Iterable[GffLine]
        An iterable of GffLine objects.

    Yields
    -------
    dict
        Dictionary corresponding to a GffLine.

    Examples
    --------
    >>> line = "TEST    test    gene    1       100     .       +       +       ID=test001;comment=Test"
    >>> list(gff2dict(read_gff(line)))  # doctest: +NORMALIZE_WHITESPACE
    [{'seqid': 'TEST', 'source': 'test', 'feature': 'gene', 'start': 1, 'end': 100, 'score': '.', 'strand': '+', 'phase': '+', 'ID': 'test001', 'comment': 'Test'}]    
    
    """

    for gff_line in gff_stream:

        d = gff_line.columns._asdict()
        d.update(gff_line.attributes)

        yield d


def gff2csv(gff_stream: Iterable[GffLine],
            file: TextIO = sys.stdout,
            write_metadata: bool = False,
            sep=',') -> None:
    
    """Converts a stream of GFFLines to a delimited stream.

    The resulting table has columns for the GFF columns 1-8 and 
    each of the unique attributes from the whole stream.
    
    Parameters
    ----------
    gff_stream : Iterable[GffLine]
        An iterable of GffLine objects.
    file : TextIO
        File-like object which has a `write` method. Default: `sys.stdout`.
    write_metadata : bool, optional
        Whether to write metadata header. Default: False.
    sep : str, optional
        Delimiter between fields. Default: ','.

    Raises
    ------
    OSError
        If the `gff_stream` is empty.

    Returns
    -------
    None

    Examples
    --------
    >>> from io import StringIO
    >>> line = '''
    ... TEST    test    gene    1       100     .       +       +       ID=test001;comment=Test
    ... TEST    test    gene    121       120     .       +       -       ID=test001;tag=test_tag
    ... '''
    >>> gff2csv(read_gff(line))  # doctest: +SKIP
    seqid,source,feature,start,end,score,strand,phase,ID,comment,tag
    TEST,test,gene,1,100,.,+,+,test001,Test,
    TEST,test,gene,121,120,.,+,-,test001,,test_tag
    >>> gff2csv(read_gff(line), sep='\\t')  # doctest: +NORMALIZE_WHITESPACE +SKIP
    seqid   source  feature start   end     score   strand  phase   ID      comment tag
    TEST    test    gene    1       100     .       +       +       test001 Test
    TEST    test    gene    121     120     .       +       -       test001         test_tag
    
    """

    gff_lines = tuple(gff_stream)
    attribute_keys = set()
    main_cols, metadata = None, None

    for i, gff_line in enumerate(gff_lines):

        if i == 0:
            main_cols = gff_line.columns._fields
            metadata = gff_line.metadata
        
        attribute_keys |= (set(gff_line.attributes))

    if main_cols is None:
        raise IOError('GFF stream is empty.')

    if write_metadata:
        
        file.write(str(metadata))

    writer = csv.DictWriter(file,
                            fieldnames=list(main_cols) + list(sorted(attribute_keys)),
                            delimiter=sep)
    writer.writeheader()
    
    for gff_line in gff2dict(gff_lines):

        writer.writerow(gff_line)

    return None