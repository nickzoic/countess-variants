""" Flexible sequence Alignment Plugin for CountESS, based on "sequence_align" """

import re
from typing import Mapping, Optional

from more_itertools import chunked

import pandas as pd

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
    description = "Flexible sequence Alignment"
    additional = "Not at all correct yet."
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
        return find_variant_string(
            self.parameters["ref"].value,
            value,
            triplet_mode = self.parameters["triplet"].value
        )
        
    def run_df(self, df: pd.DataFrame, logger: Logger) -> pd.DataFrame:

        assert isinstance(self.parameters['column'], ColumnChoiceParam)

        column = self.parameters['column'].get_column(df)

        return df.assign(variant=column.apply(self.process))
from sequence_align.pairwise import hirschberg, needleman_wunsch
from more_itertools import chunked, unzip
from itertools import groupby 

def triplets(seq):
    return [ ''.join(x) for x in chunked(seq, 3) ]

def grouper(pair):
    if not pair[0]: return 'ins'
    if not pair[1]: return 'del'
    if pair[0] == pair[1]: return 'equ'
    return 'delins'

def reversed_hirschberg(ref, seq):
    # the sequence_align library puts gaps as early in the
    # list as possible whereas we want them as 3'-ward as
    # possible.  So we just reverse the sequences, run the
    # alignment and then reverse them back.  This isn't nice
    # but it works for now.
    ref = list(reversed(ref))
    seq = list(reversed(seq))
    pairs = list(zip(*hirschberg(ref, seq, gap='', indel_score=-2)))
    pairs.reverse()
    return pairs

def variations(ref, seq):
    pairs = list(zip(*hirschberg(ref, seq, gap='', indel_score=-2)))
    #pairs = reversed_hirschberg(ref, seq)

    # adjust alignment of matches and inserts.
    # XXX this needs much more work to pass tests!
    for n in range(0, len(pairs)-1):
        if pairs[n][0] == '' and pairs[n][1] == pairs[n+1][1]:
            pairs[n] = pairs[n+1]
            pairs[n+1] = ( '', pairs[n][1] )
        elif pairs[n][1] == '' and pairs[n][0] == pairs[n+1][0]:
            pairs[n] = pairs[n+1]
            pairs[n+1] = ( pairs[n][0], '' )

    pos = 0
    for oper, group_pairs in groupby(pairs, grouper):
        rr, ss = (''.join(x) for x in unzip(group_pairs))
        pos1 = pos + 1
        pos2 = pos + len(rr)
        if oper == 'ins':
            print(f"{ref[:pos]}|{ref[pos:pos2]}|{ref[pos2:]} -> {ss}")
            if ref[pos-(pos2-pos1):pos] == ss:
                yield f"{pos1-(pos2-pos1)}_{pos1}dup"
            elif ref[pos:pos+len(ss)] == ss:
                yield f"{pos1}_{pos+len(ss)}dup"
            else:
                yield f"{pos}_{pos1}ins{ss}"
        elif oper == 'delins':
            if len(rr) == len(ss) == 1:
                yield f"{pos1}{rr}>{ss}"
            else:
                yield f"{pos1}_{pos2}delins{ss}"
        elif oper == 'del':
            if pos1 != pos2:
                yield f"{pos1}_{pos2}del"
            else:
                yield f"{pos1}del"
        pos = pos2

def find_variant_string(ref, seq, triplet_mode=False):
    
    if triplet_mode:
        vv = list(variations(triplets(ref), triplets(seq)))
    else:
        vv = list(variations(ref, seq))

    if len(vv) == 0:
        return "="
    elif len(vv) == 1:
        return vv[0]
    else:
        return '[' + ';'.join(vv) + ']'

