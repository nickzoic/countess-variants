""" Flexible sequence Alignment Plugin for CountESS, based on "sequence_align" """

import re
from typing import Mapping, Optional

from more_itertools import chunked

import pandas as pd

from sequence_align.pairwise import hirschberg

import mappy  # type: ignore

from countess.core.logger import Logger
from countess.core.parameters import (
    BooleanParam,
    ChoiceParam,
    ColumnChoiceParam,
    IntegerParam,
    StringParam,
    StringCharacterSetParam,
    FileParam,
)
from countess.core.plugins import PandasTransformPlugin

VERSION = '0.0.1'

MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


def triplets(seq):
    return [ tuple(x) for x in chunked(seq, 3) ]

class VariantsPlugin(PandasTransformPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    name = "Variant Caller"
    description = """
        Flexible sequence Alignment Plugin for CountESS, based on "sequence_align"
    """
    version = VERSION
    link = "https://github.com/CountESS-Project/countess-variants#readme"

    CHARACTER_SET = set(['A', 'C', 'G', 'T'])

    parameters = {
        "column": ColumnChoiceParam("Input Column", "sequence"),
        "ref": StringCharacterSetParam("Ref Sequence", character_set=CHARACTER_SET),
        "triplet": BooleanParam("Triplet Mode", False),
        "drop": BooleanParam("Drop Unmatched", False),
    }

    def process(self, value):
        if self.parameters["triplet"].value:
            pass
        pairs = list(zip(*hirschberg(self.parameters["seq"].value, value, gap=None)))
        
    def run_df(self, df: pd.DataFrame, logger: Logger) -> pd.DataFrame:

        assert isinstance(self.parameters['column'], ColumnChoiceParam)

        column = self.parameters['column'].get_column(df)
        prefix = self.parameters["prefix"].value

         

        return column.apply(func)
