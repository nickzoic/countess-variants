from sequence_align.pairwise import hirschberg
from more_itertools import chunked, unzip
from itertools import groupby 

def triplets(seq):
    return [ ''.join(x) for x in chunked(seq, 3) ]

def grouper(pair):
    if not pair[0]: return 'ins'
    if not pair[1]: return 'del'
    if pair[0] == pair[1]: return 'equ'
    return 'delins'

def variations(ref, seq):
    pos = 0
    pairs = zip(*hirschberg(ref, seq, gap='', indel_score=-2))
    for oper, group_pairs in groupby(pairs, grouper):
        rr, ss = (''.join(x) for x in unzip(group_pairs))
        pos1 = pos + 1
        pos2 = pos + len(rr)
        if oper == 'ins':
            if ref[pos-(pos2-pos1):pos] == ss:
                yield f"{pos1-(pos2-pos1)}_{pos1}dup"
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

