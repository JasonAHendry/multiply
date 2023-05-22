import pytest
from dataclasses import dataclass
from multiply.util.parsing import parse_parameters
from multiply.util.exceptions import DesignFileError


design_dir = "tests/fixtures/designs"

# Test exceptions
@pytest.mark.parametrize(
    "design_path, expected_exceptions",
    [
        (f"{design_dir}/doesnt-exist.ini", FileNotFoundError),
        (f"{design_dir}/wrong-suffix.txt", DesignFileError),
        (f"{design_dir}/pf-missing-comma.ini", DesignFileError),
        (f"{design_dir}/pf-missing-section.ini", DesignFileError),
        (f"{design_dir}/pf-no-genes-or-regions.ini", DesignFileError),
        (f"{design_dir}/pf-unequal-genes-ids.ini", DesignFileError),
        (f"{design_dir}/pf-invalid-genome.ini", DesignFileError)
    ],
)
def test_parse_parameters_exceptions(design_path, expected_exceptions):
    with pytest.raises(expected_exceptions):
        params = parse_parameters(design_path)


# Test non-exceptions
@dataclass
class ParsedParams:
    """
    Store expected result parsing of a design file,
    for testing
    
    """
    
    genome: str
    from_genes: bool
    from_regions: bool
    has_names: bool
    n_targets: int
    n_primer3_settings: int

@pytest.mark.parametrize(
        "design_path, result",
        [
            (f"{design_dir}/pf-default.ini", 
             ParsedParams(genome="PlasmodiumFalciparum", 
                          from_genes=True,
                          from_regions=True,
                          has_names=True,
                          n_targets=7,
                          n_primer3_settings=4)),
            (f"{design_dir}/ag-default.ini", 
             ParsedParams(genome="AnophelesGambiae", 
                          from_genes=True,
                          from_regions=False,
                          has_names=False,
                          n_targets=7,
                          n_primer3_settings=4)),

        ]      
)
def test_parse_parameters_valid(design_path, result):
    params = parse_parameters(design_path)
    assert params["genome"] == result.genome
    assert params["from_genes"] == result.from_genes
    assert params["from_regions"] == result.from_regions
    assert len(params["target_ids"]) == result.n_targets
    assert len(params["primer3_settings"]) == result.n_primer3_settings

