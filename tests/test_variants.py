from countess_variants.caller import find_variant_string

def test_single_subst():
    """
    SUBSTITUTION OF SINGLE NUCLEOTIDES (checked against Enrich2)
    """

    assert find_variant_string("AGAAGTAGAGG", "TGAAGTAGAGG") == "1A>T"
    assert find_variant_string("AGAAGTAGAGG", "AGAAGTTGTGG") == "[7A>T;9A>T]"
    assert find_variant_string("AGAAGTAGAGG", "ATAAGAAGAGG") == "[2G>T;6T>A]"


def test_duplication():
    """
    DUPLICATION OF NUCLEOTIDES

    (examples from https://varnomen.hgvs.org/recommendations/DNA/variant/duplication/ )
    """

    assert find_variant_string(
        "ATGCTTTGGTGGGAAGAAGTAGAGGACTGTTATGAAAGAGAAGATGTTCAAAAGAA",
        "ATGCTTTGGTGGGAAGAAGTTAGAGGACTGTTATGAAAGAGAAGATGTTCAAAAGAA") == "20dup"

    assert find_variant_string(
        "ATGCTTTGGTGGGAAGAAGTAGAGGACTGTTATGAAAGAGAAGATGTTCAAAAGAA",
        "ATGCTTTGGTGGGAAGAAGTAGATAGAGGACTGTTATGAAAGAGAAGATGTTCAAAAGAA") == "20_23dup"

    # (checked with Enrich2 ... it comes up with `g.6dupT` for the first
    # and g.12dupA for the second, but that trailing symbol is not how it
    # is in the HGVS docs.  The multi-nucleotide ones it gets very different
    # answers for)

    assert find_variant_string("AGAAGTAGAGG", "AGAAGTTAGAGG") == '6dup'
    assert find_variant_string("ATTGAAAAAAAATTAG", "ATTGAAAAAAAAATTAG") == '12dup'
    assert find_variant_string("AGAAGTAGAGG", "AGAAGTAGATAGAGG") == '6_9dup'
    assert find_variant_string("AAAACTAAAA", "AAAACTCTAAAA") == '5_6dup'
    assert find_variant_string("AAAACTGAAAA", "AAAACTGCTGAAAA") == '5_7dup'
    assert find_variant_string("AAAACTGAAAA", "AAAACTGTGAAAA") == '6_7dup'

    # It should always pick the 3'-most copy to identify as
    # duplicate

    assert find_variant_string("ATGGCCGCCAGCCAA", "ATGGCCGCCGCCAGCCAA") == '7_9dup'

def test_insertion():
    """
    INSERTION OF NUCLEOTIDES
    """

    assert find_variant_string("AGAAGTAGAGG", "AGAAGTCAGAGG") == '6_7insC'
    assert find_variant_string("AGAAGTAGAGG", "AGAAGTCATAGAGG") == '6_7insCAT'

    #        0        1
    #        12345678901234
    #        AGAAGTAGAGG
    #        AGAAGTGAAAGAGG
    #         ^^^  ^^^
    #
    # GAA appears earlier in the sequence, but it's too short to bother
    # including as a reference.

    assert find_variant_string("AGAAGTAGAGG", "AGAAGTGAAAGAGG") == '6_7insGAA'

    # the sequence CTTTTTTTTT is long enough to use a reference for.

    assert find_variant_string("ACTTTTTTTTTAA", "ACTTTTTTTTTACTTTTTTTTTA") == '12_13ins2_11'

    # the most 3'ward copy of the sequence should be the one which is
    # referred to.

    assert find_variant_string(
        "GGCATCATCATCATGGCATCATCATCATGGGG", 
        "GGCATCATCATCATGGCATCATCATCATGGCATCATCATCATGG") == '30_31ins17_28'

    assert find_variant_string(
        "GGCATCATCATCATGGCATCATCATCATGGGGCATCATCATCATGG",
        "GGCATCATCATCATGGCATCATCATCATGGCATCATCATCATGGCATCATCATCATGG") == '30_31ins33_44'

    # (example from Q&A on https://varnomen.hgvs.org/recommendations/DNA/variant/insertion/ )
    # "How should I describe the change ATCGATCGATCGATCGAGGGTCCC to
    # ATCGATCGATCGATCGAATCGATCGATCGGGTCCC? The fact that the inserted sequence
    # (ATCGATCGATCG) is present in the original sequence suggests it derives
    # from a duplicative event."
    # "The variant should be described as an insertion; g.17_18ins5_16."
    #
    #         0        1         2         3
    #         1234567890123456789012345678901234
    #         ATCGATCGATCGATCGAGGGTCCC
    #         ATCGATCGATCGATCGAATCGATCGATCGGGTCCC
    #             ^^^^^^^^^^^  ^^^^^^^^^^^
    #
    # right so the duplicated part is 5_15 and its inserted between 17 and 18.
    # I think the example in the Q&A is wrong.

    assert find_variant_string("ATCGATCGATCGATCGAGGGTCCC", "ATCGATCGATCGATCGAATCGATCGATCGGGTCCC") == '17_18ins5_15'

def test_deletion():
    """
    DELETION OF NUCLEOTIDES (checked against Enrich2)
    enrich2 has g.5_5del for the first one which isn't correct though
    """

    assert find_variant_string("AGAAGTAGAGG", "AGAATAGAGG") == '5del'
    assert find_variant_string("AGAAGTAGAGG", "AGAAAGAGG") == '5_6del'

    # MULTIPLE VARIATIONS

    assert find_variant_string("GATTACA", "GATTACA") == '='
    assert find_variant_string("GATTACA", "GTTTACA") == '2A>T'
    assert find_variant_string("GATTACA", "GTTTAGA") == '[2A>T;6C>G]'
    assert find_variant_string("GATTACA", "GTTCAGA") == '[2A>T;4T>C;6C>G]'

    assert find_variant_string("ATGGTTGGTTC", "ATGGTTGGTGGTTC") == '7_9dup'
    assert find_variant_string("ATGGTTGGTTCG", "ATGGTTGGTGGTTC") == '[7_9dup;12del]'
    assert find_variant_string("ATGGTTGGTTC", "ATGGTTGGTGGTTCG") == '[7_9dup;11_12insG]'
