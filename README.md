# 💻 bioino

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/bioino/python-publish.yml)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bioino)
![PyPI](https://img.shields.io/pypi/v/bioino)

Command-line tools and Python API for interconverting FASTA, GFF, and CSV.

**bioino** converts tables to FASTA, and GFF to tables. It also provides a Python API
for reading, writing, and querying GFF and FASTA files.

## Installation

### The easy way

```bash
pip install bioino
```

### From source

Clone the repository, then `cd` into it and run:

```bash
pip install -e .
```

## Usage

### Command line

Info goes to `stderr`, so output can be piped freely.

#### `gff2table`

Convert a GFF file to TSV (default) or CSV.

```bash
$ printf 'test_seq\ttest_source\tgene\t1\t10\t.\t+\t.\tID=test01;attr1=+\n' \
    | bioino gff2table 2>/dev/null
seqid   source  feature start   end     score   strand  phase   ID      attr1
test_seq        test_source     gene    1       10      .       +       .       test01  +
```

```bash
$ printf 'test_seq\ttest_source\tgene\t1\t10\t.\t+\t.\tID=test01;attr1=+\n' \
    | bioino gff2table -f CSV 2>/dev/null
seqid,source,feature,start,end,score,strand,phase,ID,attr1
test_seq,test_source,gene,1,10,.,+,.,test01,+
```

Pass `--metadata` / `-m` to include the GFF header as commented lines in the output.

#### `table2fasta`

Convert a CSV or TSV table of sequences to FASTA.

```bash
$ printf 'name\tseq\tdata\nSeq1\tAAAAA\tSome-info\n' \
    | bioino table2fasta -n name -s seq -d data 2>/dev/null
>Seq1 data=Some-info
AAAAA
```

Multiple `--name` columns are concatenated with `_`; multiple `--description` columns
are formatted as `key=value` pairs separated by `;`.

#### Detailed usage

```
usage: bioino [-h] [--version] {gff2table,table2fasta} ...

Interconvert some bioinformatics file formats.

options:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit

Sub-commands:
  {gff2table,table2fasta}
    gff2table           Convert a GFF to a TSV file.
    table2fasta         Convert a CSV or TSV of sequences to a FASTA file.
```

```
usage: bioino gff2table [-h] [--format {TSV,CSV}] [--metadata]
                        [--output OUTPUT]
                        [input]

positional arguments:
  input                 Input file in GFF format. Default: stdin.

options:
  -h, --help            show this help message and exit
  --format {TSV,CSV}, -f {TSV,CSV}
                        Output format. Default: "TSV".
  --metadata, -m        Write GFF header as commented lines.
  --output OUTPUT, -o OUTPUT
                        Output file. Default: stdout.
```

```
usage: bioino table2fasta [-h] [--format {TSV,CSV}] [--sequence SEQUENCE]
                          --name [NAME ...] [--description [DESCRIPTION ...]]
                          [--worksheet WORKSHEET] [--output OUTPUT]
                          [input]

positional arguments:
  input                 Input table file (TSV, CSV, or XLSX). Default: stdin.

options:
  -h, --help            show this help message and exit
  --format {TSV,CSV}, -f {TSV,CSV}
                        Input format. Default: "TSV".
  --sequence SEQUENCE, -s SEQUENCE
                        Column to take sequence from. Default: "sequence".
  --name [NAME ...], -n [NAME ...]
                        Column(s) for sequence name. Concatenated with "_",
                        spaces replaced with "-". Required.
  --description [DESCRIPTION ...], -d [DESCRIPTION ...]
                        Column(s) for sequence description. Formatted as
                        "key=value" pairs separated by ";", spaces replaced
                        with "_". Default: omitted.
  --worksheet WORKSHEET, -w WORKSHEET
                        For XLSX files, the worksheet to read. Default: "Sheet 1".
  --output OUTPUT, -o OUTPUT
                        Output file. Default: stdout.
```

---

### Python API

#### FASTA

`FastaSequence` is a dataclass holding a sequence name, description, and sequence string.
`FastaCollection` wraps an iterable of `FastaSequence` objects.

```python
>>> from bioino import FastaSequence, FastaCollection

>>> seq1 = FastaSequence("example", "This is a description", "ATCG")
>>> seq2 = FastaSequence("example2", "This is another sequence", "GGGAAAA")
>>> FastaCollection([seq1, seq2]).write()
>example This is a description
ATCG
>example2 This is another sequence
GGGAAAA
```

Read from a file handle or filename with `FastaCollection.from_file()`:

```python
>>> from io import StringIO
>>> buf = StringIO()
>>> FastaCollection([seq1, seq2]).write(buf)
>>> buf.seek(0)
0
>>> FastaCollection.from_file(buf).write()
>example This is a description
ATCG
>example2 This is another sequence
GGGAAAA
```

Build a `FastaCollection` from a Pandas DataFrame with `FastaCollection.from_pandas()`.
The `names` columns are concatenated with `name_sep` (default `_`); `descriptions` columns
are formatted as `key=value` pairs separated by `desc_sep` (default `;`).

```python
>>> import pandas as pd
>>> from bioino import FastaCollection

>>> df = pd.DataFrame(dict(
...     seq=['atcg', 'aaaa'],
...     title=['seq1', 'seq2'],
...     info=['SeqA', 'SeqB'],
...     score=[1, 2],
... ))
>>> FastaCollection.from_pandas(df, sequence='seq',
...                             names=['title'],
...                             descriptions=['info', 'score']).write()
>seq1 info=SeqA;score=1
atcg
>seq2 info=SeqB;score=2
aaaa
>>> FastaCollection.from_pandas(df, sequence='seq',
...                             names=['title', 'info'],
...                             descriptions=['score']).write()
>seq1_SeqA score=1
atcg
>seq2_SeqB score=2
aaaa
```

---

#### GFF

Makes an attempt to conform to GFF3 but makes no guarantees.

**Reading and writing**

`GffFile.from_file()` accepts a file handle or filename and returns a `GffFile` that
streams records lazily.

```python
>>> from io import StringIO
>>> from bioino import GffFile

>>> lines = [
...     "##meta1 item1",
...     "#meta2  item2  comment",
...     "\t".join("test_seq test_source gene 1 10 . + . ID=test01;attr1=+".split()),
...     "\t".join("test_seq test_source gene 9 100 . + . Parent=test01;attr2=+".split()),
... ]
>>> gff = GffFile.from_file(StringIO("\n".join(lines)))
>>> gff.write()
##meta1 item1
#meta2  item2  comment
test_seq    test_source     gene    1       10      .       +       .       ID=test01;attr1=+
test_seq    test_source     gene    9       100     .       +       .       Parent=test01;attr2=+
```

**Converting to table**

`GffFile.to_csv()` writes a flat table with one row per GFF line, columns for the eight
standard GFF fields plus all unique attribute keys. Use `sep='\t'` for TSV output.

```python
>>> from io import StringIO
>>> from bioino import GffFile

>>> lines = [
...     "\t".join("TEST test gene 1 100 . + + ID=test001;comment=Test".split()),
...     "\t".join("TEST test gene 121 120 . + - ID=test001;tag=test_tag".split()),
... ]
>>> GffFile.from_file(StringIO("\n".join(lines))).to_csv()
seqid,source,feature,start,end,score,strand,phase,ID,comment,tag
TEST,test,gene,1,100,.,+,+,test001,Test,
TEST,test,gene,121,120,.,+,-,test001,,test_tag
```

**Interconversion**

`GffLine.from_dict()` constructs a `GffLine` from a dictionary. Keys matching the standard
GFF column names (`seqid`, `source`, `feature`, `start`, `end`, `score`, `strand`, `phase`)
populate the columns; all other keys become attributes.

```python
>>> from bioino import GffLine

>>> d = dict(seqid='TEST', source='test', feature='gene',
...          start=1, end=100, score='.', strand='+', phase='+')
>>> print(GffLine.from_dict(d))
TEST    test    gene    1       100     .       +       +

>>> d.update(dict(ID='test001', comment='This is a test'))
>>> GffLine.from_dict(d).write()
TEST    test    gene    1       100     .       +       +       ID=test001;comment=This is a test
```

`GffFile.as_dict()` yields each line as a flat dictionary:

```python
>>> from io import StringIO
>>> from bioino import GffFile

>>> lines = [
...     "TEST\ttest\tgene\t1\t100\t.\t+\t+\tID=test001;comment=Test",
...     "TEST2\ttest2\tgene\t101\t200\t.\t+\t+\tID=test002;comment=Test2",
... ]
>>> list(GffFile.from_file(StringIO("\n".join(lines))).as_dict())
[{'seqid': 'TEST', 'source': 'test', 'feature': 'gene', 'start': 1, 'end': 100,
  'score': '.', 'strand': '+', 'phase': '+', 'ID': 'test001', 'comment': 'Test'},
 {'seqid': 'TEST2', 'source': 'test2', 'feature': 'gene', 'start': 101, 'end': 200,
  'score': '.', 'strand': '+', 'phase': '+', 'ID': 'test002', 'comment': 'Test2'}]
```

**Positional lookup**

`GffFile` can build a per-chromosome interval index for fast positional annotation queries.
Pass `lookup=True` to `GffFile.from_file()`.

```python
>>> from io import StringIO
>>> from bioino import GffFile

>>> lines = [
...     "\t".join(["chr1", "src", "gene", "10",  "50",  ".", "+", ".", "ID=g1;Name=geneA"]),
...     "\t".join(["chr1", "src", "gene", "100", "150", ".", "+", ".", "ID=g2;Name=geneB"]),
...     "\t".join(["chr2", "src", "gene", "20",  "80",  ".", "-", ".", "ID=g3;Name=geneC"]),
... ]
>>> gff = GffFile.from_file(StringIO("\n".join(lines)), lookup=True)
```

Query with `lookup_at(seqid, pos)`, which returns a tuple of `GffLine` objects covering
that position. Each returned line has `locus_tag` and `offset` attributes computed for
that exact position.

```python
# Gene body — offset from annotated start (+ strand) or end (- strand)
>>> r = gff.lookup_at('chr1', 30)
>>> r[0].attributes['locus_tag'], r[0].attributes['offset']
('geneA', 20)

# Intergenic — first half of gap attributed to upstream gene
>>> r = gff.lookup_at('chr1', 75)
>>> r[0].attributes['locus_tag'], r[0].attributes['offset']
('_down-geneA', 65)

# Intergenic — second half of gap attributed to downstream gene
>>> r = gff.lookup_at('chr1', 76)
>>> r[0].attributes['locus_tag'], r[0].attributes['offset']
('_up-geneB', 24)

# Up to 1000 bp past the last annotated feature is covered
>>> r = gff.lookup_at('chr1', 200)
>>> r[0].attributes['locus_tag'], r[0].attributes['offset']
('_down-geneB', 100)

# Each chromosome is indexed independently
>>> r = gff.lookup_at('chr2', 50)
>>> r[0].attributes['locus_tag'], r[0].attributes['offset']
('geneC', 30)

# Returns an empty tuple for unknown seqids or positions outside all intervals
>>> gff.lookup_at('chrX', 50)
()
```

The lookup index:
- handles multi-chromosome GFFs
- only indexes parent features (`Name` attribute present, no `Parent` attribute)
- ignores feature types `region` and `repeat_region`
- stores references to the original `GffLine` objects; offsets are computed on demand

---

## Suggestions, issues, fixes

File an issue [here](https://github.com/scbirlab/bioino/issues).

## Documentation

API reference at [bioino.readthedocs.org](https://bioino.readthedocs.org).
