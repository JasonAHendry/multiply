import pytest
from multiply.util.parsing import parse_parameters


# TODO:
# - This is pretty useless testing
# - We would want to check that the *expected* exceptions are thrown
design_dir = "tests/fixtures/designs"
@pytest.mark.parametrize(
    "design_path, expected_exceptions",
    [
        (f"{design_dir}/doesnt-exist.ini", FileNotFoundError),
        (f"{design_dir}/wrong-suffix.txt", ValueError),
        (f"{design_dir}/pf-missing-comma.ini", ValueError),
        (f"{design_dir}/pf-missing-section.ini", ValueError),
        (f"{design_dir}/pf-no-genes-or-regions.ini", ValueError),
        (f"{design_dir}/pf-unequal-genes-ids.ini", ValueError)
    ],
)
def test_parse_parameters_exceptions(design_path, expected_exceptions):
    with pytest.raises(expected_exceptions):
        params = parse_parameters(design_path)

