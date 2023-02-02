import pytest
import math
from itertools import combinations, product
from multiply.align.algorithms import PrimerDimerLike


# PREPARE THE ALIGNER
aligner = PrimerDimerLike()
aligner.load_parameters()


# TESTS
combos = [(b, l) for b, l in product("ATCG", [5, 6, 7, 8, 9])]  # can test everything > `end_length`
@pytest.mark.parametrize(
    "base, length", combos
)
def test_homopolymer_primer_with_overlap_scores(base, length):
    """
    Test the dimer scores for homopolymers;
    this is pretty easy because the best alignment is easy
    to predict
    
    """

    # Get bases
    rc_base = aligner.rc_map[base]

    # Create primers
    primer1 = base * length
    primer2 = rc_base * length
    aligner.set_primers(
        primer1=primer1,
        primer2=primer2,
        primer1_name="primer1",
        primer2_name="primer2"
    )

    # Compute expected score
    # Should be the longest match allowing 5' overhangs on
    # both sides.
    match_score = aligner.nn_scores[f"{base*2}/{rc_base*2}"]
    matches_score = match_score * (length - 2)
    end_scores = 2 * aligner.end_length * aligner.end_bonus
    total_score = matches_score + end_scores

    # Align
    aligner.align()

    print(aligner.get_alignment_string())

    # Check we have the expected score
    assert math.isclose(aligner.score, total_score)


combos = list(product(combinations("ACG", 2), [4, 5, 6, 7]))  # AT will mismatch to produce overhangs
@pytest.mark.parametrize(
    "bases, n", combos
)
def test_two_homopolymers_without_overlap_scores(bases, n):
    """
    The following cases...

    5-AAAACCCC-3
    3-TTTTGGGG-5

    ...are easy to calculate.
    
    """

    # Create primers
    base1, base2 = bases
    rbase1 = aligner.rc_map[base1]
    rbase2 = aligner.rc_map[base2]
    primer1 = base1*n + base2*n
    primer2 = rbase2*n + rbase1*n

    # Run alignment
    aligner.set_primers(
        primer1=primer1,
        primer2=primer2,
        primer1_name="test1",
        primer2_name="test2"
    )
    aligner.align()

    # Compute expected scores
    a = aligner.nn_scores[f"{base1*2}/{rbase1*2}"]*(n-1)
    b = aligner.nn_scores[f"{base1+base2}/{rbase1+rbase2}"]
    c = aligner.nn_scores[f"{base2*2}/{rbase2*2}"]*(n-1)
    total_score = a + b + c

    # Check
    assert math.isclose(aligner.score, total_score)


# def test_no_end_score():
#     """
#     """
#     pass


# def test_end_scoring():
#     """
#     Test the dimer scores for homopolymers;
#     this is pretty easy because the best alignment is easy
#     to predict
    
#     """
#     pass
