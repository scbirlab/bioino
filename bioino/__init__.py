from importlib.metadata import version

appname = "bioino"
__version__ = version(appname)

from .fasta import *
from .gff.column import GffColumns
from .gff.gff import GffFile
from .gff.line import GffLine
