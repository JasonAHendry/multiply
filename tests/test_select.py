import pytest
from multiply.select.multiplex import Multiplex

# Fixtures
m1 = Multiplex(cost=-1, primer_pairs=["A", "C", "B"])
m2 = Multiplex(cost=-2, primer_pairs=["D", "F", "E"])
m3 = Multiplex(cost=-3, primer_pairs=["A", "B", "D"])
m4 = Multiplex(cost=8, primer_pairs=["A", "C", "B"]) # same pairs as m1, but different cost
ms = [m1, m2, m3, m1, m1, m2, m4]

def test_multiplex_primer_sort():
    """
    Test that the class internally sorts
    primer names
    """
    assert m1.primer_pairs == ("A", "B", "C")
    assert m2.primer_pairs == ("D", "E", "F")

def test_multiplex_unique():
    """
    Test that the class correctly reduces
    to a unique set of multiplexes

    """
    uniq_ms = [m1, m2, m3]
    assert len(set(ms)) == len(uniq_ms)
    assert set(ms) == set(uniq_ms)

def test_multiplex_sort_and_unique():
    ms_sorted = sorted(set(ms))
    assert ms_sorted == [m1, m3, m2]

    