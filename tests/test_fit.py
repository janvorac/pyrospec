import json
import numpy as np
import pytest
from pyrospec import fit, Params, estimate_intensity, SpectralSystemEnum


@pytest.mark.parametrize(
    "input_file",
    [
        "tests/inputs/fit_input_x_shifted.json",
    ],
)
def test_fit_with_sample_data(input_file):
    """Test the fit function with sample data from different input files."""
    # Load test data
    with open(input_file, "r") as f:
        input_data = json.load(f)

    y_measured = np.array(input_data["y"])
    initial_guess = Params(**input_data["initialGuess"])

    ohax_estimated_intensity = estimate_intensity(
        initial_guess=initial_guess,
        species=SpectralSystemEnum.OHAX,
        y_measured=y_measured,
    )

    initial_guess.sim.species_params[SpectralSystemEnum.OHAX].intensity = (
        ohax_estimated_intensity
    )

    # Run the fit
    params, res = fit(
        y_measured=y_measured,
        initial_guess=initial_guess,
        settings=None,
    )

    # Basic assertions to ensure the fit ran successfully
    assert params is not None
    assert res is not None
    assert hasattr(
        res, "success"
    )  # scipy optimization result should have success attribute
    assert res.success is True
