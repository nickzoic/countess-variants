from countess_variants.caller import find_variant_string

def test_triplet_subst():
    """
    SUBSTITUTION OF SINGLE NUCLEOTIDES (checked against Enrich2)
    """

    assert find_variant_string("ATGGTT", "ATGGGTGTT", triplet_mode=True) == "3_4insGGT"


