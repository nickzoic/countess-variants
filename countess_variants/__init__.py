""" Flexible sequence Alignment Plugin for CountESS, based on "sequence_align" """

import re
from typing import Mapping, Optional

import pandas as pd
from countess.core.logger import Logger
from countess.core.parameters import (BooleanParam, ColumnChoiceParam,
                                      StringCharacterSetParam)
from countess.core.plugins import PandasTransformPlugin
from more_itertools import chunked

from .caller import find_variant_string

VERSION = "0.0.1"


class VariantsPlugin(PandasTransformPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    name = "Variant Caller"
    description = "Flexible sequence Alignment"
    additional = "Not at all correct yet."
    version = VERSION
    link = "https://github.com/CountESS-Project/countess-variants#readme"

    CHARACTER_SET = set(["A", "C", "G", "T"])

    parameters = {
        "column": ColumnChoiceParam("Input Column", "sequence"),
        "ref": StringCharacterSetParam("Ref Sequence", character_set=CHARACTER_SET),
        "triplet": BooleanParam("Triplet Mode", False),
        "drop": BooleanParam("Drop Unmatched", False),
    }

    def process(self, value):
        """process a single sequence into a variant string"""
        return find_variant_string(
            self.parameters["ref"].value, value, triplet_mode=self.parameters["triplet"].value
        )

    def run_df(self, df: pd.DataFrame, logger: Logger) -> pd.DataFrame:
        """process dataframe"""
        assert isinstance(self.parameters["column"], ColumnChoiceParam)

        column = self.parameters["column"].get_column(df)

        return df.assign(variant=column.apply(self.process))
