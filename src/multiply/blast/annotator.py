from collections import namedtuple


# ================================================================================
# Annotations to add to BLAST results
# - Always add from_3prime
# - Two additional metrics used: evalue and alignment length
# - Thresholds should be available in /settings/blast
# - Changing this dictionary changes all downstream behaviour, no code modification needed
# ================================================================================


ANNOTATIONS = {
    "from_3prime": lambda row: row["qend"] == row["length"],
    "over12_alignment_3prime": lambda row: row["length"] >= 12
    and row["pident"] == 100
    and row["from_3prime"],
    "below4_evalue_3prime": lambda row: row["evalue"] < 4 and row["from_3prime"],
    "predicted_bound": lambda row: row["over12_alignment_3prime"]
    or row["below4_evalue_3prime"],
}


# ================================================================================
# BLAST output information for a primer pair
#
#
# ================================================================================


PrimerBlastRecord = namedtuple(
    "PrimerBlastRecord",
    ["primer_name", "primer_pair_name", "target_name", "total_alignments"]
    + list(ANNOTATIONS),
)
