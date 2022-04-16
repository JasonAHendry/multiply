import pytest
from multiply.util.parsing import parse_parameters
from multiply.util.exceptions import DesignFileError


# TODO:
# - Have tested the exceptions are thrown
# - Should probably test that params loaded correctly for non-exception cases

design_dir = "tests/fixtures/designs"
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

