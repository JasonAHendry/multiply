import numpy as np
from collections import Counter


def get_homopolymer_runs(seq, l_max=None):
    """
    For a given sequence `seq`, produce a homopolymer
    block size encoding

    i.e.
    ATTCCC = 1, 2, 2, 3, 3, 3

    """

    # Store
    h = np.zeros(len(seq), "int8")
    h[:] = -1

    # Initialise
    i = 0
    p = seq[i]
    l = 0

    # Iterate
    for c in seq:
        if c == p:
            l += 1
        else:
            h[i : (i + l)] = l
            i += l
            p = c
            l = 1

    # Terminate
    h[i : (i + l)] = l
    if l_max is not None:
        h[h > l_max] = l_max
        assert (h <= l_max).all(), "Error in homopolymer encoding."
    assert (h > 0).all(), "Error in homopolymer encoding."

    return h


def calc_sliding_percentGC(seq, window):
    """
    Calculate GC content in a sliding
    window over the sequence
    
    """
    # Preapre
    n = len(seq)
    gc = np.zeros(n)
    
    # Iterate
    for i in range(n - window + 1):
        sub = seq[i:(i+window)]
        cnts = Counter(sub)
        gc[i] = cnts["G"] + cnts["C"]
        
    gc /= window
        
    return gc


def get_array_encoding(seq):
    """
    Convert `seq` into an integer array

    """

    n = len(seq)
    a = np.zeros((4, n))
    dt = {"A" : 0, "T" : 1, "C" : 2, "G" : 3, "[":np.nan, "]":np.nan}
    for i, base in enumerate(seq):
        e = dt[base]
        if e == e:
            a[e, i] = 1

    return a