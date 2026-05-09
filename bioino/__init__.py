from importlib.metadata import version

appname = "bioino"
__version__ = version(appname)

from .fasta import *
from .gff import *