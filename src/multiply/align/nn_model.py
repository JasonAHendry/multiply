import json
from itertools import product


def create_nn_score_dt(match_json, single_mismatch_json, double_mismatch_score=0.2):
    """
    Create the Gibb's free energy nearest neighbour scoring dictionary

    """
    # Load match and single mismatch .jsons
    match_dt = json.load(open(match_json, "r"))
    single_mismatch_dt = json.load(open(single_mismatch_json, "r"))

    # Set all as double mismatches; then update
    nts = ["A", "T", "C", "G"]
    nn_score_dt = {
        "".join(watson) + "/" + "".join(crick): double_mismatch_score
        for watson in product(nts, repeat=2)
        for crick in product(nts, repeat=2)
    }

    # Update
    nn_score_dt.update(match_dt)
    nn_score_dt.update(single_mismatch_dt)

    return nn_score_dt
