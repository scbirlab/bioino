"""Utilities for working with GFF files."""
from typing import Dict, Iterable, Optional, Tuple, Union
from io import TextIOWrapper

from carabiner.cast import cast


def _cast_to_file_handle(f: Union[str, TextIOWrapper]) -> TextIOWrapper:
    if isinstance(f, TextIOWrapper):
        return f
    else:
        return cast(f, to=TextIOWrapper)
