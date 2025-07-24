import numpy as np

from models import MeasuredSpectrumParams, SimulatedSpectrumParams, SpectralSystemParams
from enums import SpectralSystemEnum
from pyrospec import spectrum
from specdata import generate_spectrum, SpecDB
from scipy.optimize import minimize


def get_residuals(
    measured: np.ndarray,
    measured_spectrum_params: MeasuredSpectrumParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    simulations: dict[SpectralSystemEnum, SpecDB],
) -> np.ndarray:

    x_measured = _get_x_measured(measured_spectrum_params, len(measured))

    simulated_spec: spectrum.Spectrum = generate_spectrum(
        simspec_params=simulated_spectrum_params,
        step=measured_spectrum_params.wav_step,
        wmin=measured_spectrum_params.wav_start,
        wmax=np.max(x_measured),
        simulations=simulations,
    )

    measured_spec = spectrum.Spectrum(x=x_measured, y=measured)

    sim_spec_matched, exp_spec = spectrum.match_spectra(
        sim_spec=simulated_spec, exp_spec=measured_spec
    )

    return sim_spec_matched.y - exp_spec.y


def _get_x_measured(
    measured_spectrum_params: MeasuredSpectrumParams, length: int
) -> np.ndarray:
    return (
        measured_spectrum_params.wav_start
        + np.arange(length) * measured_spectrum_params.wav_step
        + measured_spectrum_params.wav_2nd * np.arange(length) ** 2
    )


def _min_func(
    x0,
    y_measured,
    measured_spectrum_params: MeasuredSpectrumParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    x0_meaning,
):
    _update_params(measured_spectrum_params, simulated_spectrum_params, x0, x0_meaning)
    residuals = get_residuals(
        measured=y_measured,
        measured_spectrum_params=measured_spectrum_params,
        simulated_spectrum_params=simulated_spectrum_params,
        simulations={
            # SpectralSystemEnum.OHAX: SpecDB("OHAX.db"),
            system: SpecDB(system.filename)
            for system in SpectralSystemEnum  # TODO: think if we really need all systems always
        },
    )
    # return residuals
    return np.sum(residuals**2)  # return sum of squares for minimization


def _update_params(
    measured_spectrum_params: MeasuredSpectrumParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    x0,
    x0_meaning,
):
    species_param_keywords = ["Tvib", "Trot", "intensity"]
    for param, value in zip(x0_meaning, x0):
        if param in measured_spectrum_params.__dict__:
            measured_spectrum_params.__dict__[param] = value
        elif param in simulated_spectrum_params.__dict__:
            simulated_spectrum_params.__dict__[param] = value
        elif any(kw in param for kw in species_param_keywords):
            species_name = param.split("|")[0]
            if species_name in simulated_spectrum_params.species_params:
                species_param = param.split("|")[1]
                simulated_spectrum_params.species_params[species_name].__dict__[
                    species_param
                ] = value
            else:
                raise ValueError(f"Species {species_name} not found in params.")


if __name__ == "__main__":
    import json

    with open("/home/jan/packages/pyrospec/tests/inputs/simulate_input.json", "r") as f:
        input_data = json.load(f)
    y_measured = np.array(input_data["y"])

    initial_measured_params = MeasuredSpectrumParams(
        wav_start=305.55,
        wav_step=0.0318,
        wav_2nd=-2e-7,
    )

    initial_simulated_params = SimulatedSpectrumParams(
        slitf_gauss=0.0396,
        slitf_lorentz=0.0138,
        baseline=0.55,
        baseline_slope=0,
        species_params={
            SpectralSystemEnum.OHAX: SpectralSystemParams(
                Tvib=3000,
                Trot=3000,
                intensity=1,
            )
        },
    )

    x0 = [3000, 3000, 1]
    x0_meaning = ["OHAX|Tvib", "OHAX|Trot", "OHAX|intensity"]

    res = minimize(
        fun=_min_func,
        x0=x0,
        args=(
            y_measured,
            initial_measured_params,
            initial_simulated_params,
            x0_meaning,
        ),
        method="Nelder-Mead",
        # options={'disp': True}
    )

    print(res)
