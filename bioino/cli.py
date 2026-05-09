"""Command-line interface to bioino."""

from typing import Callable

from argparse import FileType, Namespace
import sys

from carabiner.cliutils import CLIApp, CLICommand, CLIOption, clicommand
from carabiner.pd import read_table

from .fasta import FastaCollection
from .gff.gff import GffFile

from . import appname, __version__


def _allow_broken_pipe(f: Callable) -> Callable:
    def _f(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except BrokenPipeError:
            pass
    return _f


@clicommand(message='Converting GFF to table with the following parameters')
def _gff2table(args: Namespace) -> None:

    gff_file = GffFile.from_file(args.input)
    _allow_broken_pipe(gff_file.to_csv)(
        args.output,
        write_metadata=args.metadata,
        sep={"tsv": "\t", "csv": ","}[args.format.casefold()],
    )
    return None


@clicommand(message='Generating FASTA from tables with the following parameters')
def _table2fasta(args: Namespace) -> None:
    table = read_table(
        args.input, 
        format=args.format,
    )
    
    fasta_collection = FastaCollection.from_file(
        args.input,
    )
    _allow_broken_pipe(fasta_collection.write)(file=args.output)
    return None


def main() -> None:

    inputs = CLIOption('input', 
                        type=FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='Input file in GFF format.')
    format = CLIOption('--format', '-f', 
                        type=str,
                        default='TSV',
                        choices=['TSV', 'CSV'],
                        help='File format.')
    metadata = CLIOption('--metadata', '-m', 
                        action='store_true',
                        help='Write GFF header as commented lines.')
    outputs = CLIOption('--output', '-o', 
                        type=FileType('w'),
                        default=sys.stdout,
                        help='Output file. Default: STDOUT')
    
    sequence = CLIOption('--sequence', '-s', 
                        type=str, 
                        default='sequence',
                        help='Column to take sequence from.')
    name = CLIOption('--name', '-n', 
                        nargs='*',
                        type=str,
                        required=True,
                        help='Column(s) to take sequence name from. '
                             'Concatenates values with "_", '
                             'replaces spaces with "-".')
    description = CLIOption('--description', '-d', 
                        nargs='*',
                        type=str, 
                        default=[],
                        help='Column(s) to take sequence description from. '
                             'Concatenates values with ";", '
                             'replaces spaces with "_". Default: don\'t use.')
    worksheet = CLIOption('--worksheet', '-w', 
                          type=str, default='Sheet 1',
                          help='For XLSX files, the worksheet to take the table from.')

    gff2table = CLICommand("gff2table",
                           description="Convert a GFF to a TSV file.",
                           main=_gff2table,
                           options=[inputs, format, metadata, outputs])
    table2fasta = CLICommand("table2fasta",
                             description="Convert a CSV or TSV of sequences to a FASTA file.",
                             main=_table2fasta,
                             options=[inputs, format, sequence, name, description, worksheet, outputs])

    app = CLIApp(
        appname,
        version=__version__,
        description="Interconvert some bioinformatics file formats.",
        commands=[gff2table, table2fasta],
    )
    
    app.run()

    return None


if __name__ == '__main__':

    main()