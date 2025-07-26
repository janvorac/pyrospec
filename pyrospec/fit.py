import copy

import numpy as np
from .enums import SpectralSystemEnum
from .models import (
    FitSettings,
    MeasuredSpectrumParamFitSettings,
    MeasuredSpectrumParams,
    Params,
    SimulatedSpectrumParamFitSettings,
    SimulatedSpectrumParams,
    SpectralSystemParamFitSettings,
    SpectralSystemParams,
)
from scipy.optimize import Bounds, minimize
from .specdata import SpecDB, generate_spectrum

from . import spectrum


def fit(
    y_measured: np.ndarray, initial_guess: Params, settings: FitSettings | None = None
):
    """
    Fit the measured spectrum to the simulated spectrum using the provided initial guess and settings.

    :param y_measured: Measured spectrum data.
    :param initial_guess: Initial parameters for the fit.
    :param settings: Optional fit settings to control parameter bounds and fixed states.
    :return: Fitted parameters and optimization result.
    """
    params = copy.deepcopy(initial_guess)

    measured_spectrum_params = params.meas
    simulated_spectrum_params = params.sim

    x0 = []
    x0_meaning = []

    if settings is None:
        settings = _make_default_fit_settings(initial_guess)

    bounds = []

    for param_name, param_value in measured_spectrum_params.__dict__.items():
        if hasattr(settings.meas, f"{param_name}_fixed") and getattr(
            settings.meas, f"{param_name}_fixed"
        ):
            # Only not fixed parameters go into the optimization
            continue
        x0.append(param_value)
        x0_meaning.append(param_name)

        lb = getattr(settings.meas, f"{param_name}_min", -np.inf)
        ub = getattr(settings.meas, f"{param_name}_max", np.inf)
        bounds.append((lb, ub))

    for param_name, param_value in simulated_spectrum_params.__dict__.items():
        if (
            hasattr(settings.sim, f"{param_name}_fixed")
            and getattr(settings.sim, f"{param_name}_fixed")
        ) or param_name == "species_params":
            continue
        x0.append(param_value)
        x0_meaning.append(param_name)
        lb = getattr(settings.sim, f"{param_name}_min", -np.inf)
        ub = getattr(settings.sim, f"{param_name}_max", np.inf)
        bounds.append((lb, ub))

    # Add species parameters if available
    if simulated_spectrum_params.species_params:
        for species, species_params in simulated_spectrum_params.species_params.items():
            for param_name, param_value in species_params.__dict__.items():
                if hasattr(
                    settings.species[species], f"{param_name}_fixed"
                ) and getattr(settings.species[species], f"{param_name}_fixed"):
                    continue
                x0.append(param_value)
                x0_meaning.append(f"{species.name}|{param_name}")
                lb = getattr(
                    settings.species[species],
                    f"{param_name}_min",
                    -np.inf,
                )
                ub = getattr(
                    settings.species[species],
                    f"{param_name}_max",
                    np.inf,
                )
                bounds.append((lb, ub))

    # Perform minimization
    res = minimize(
        fun=_min_func,
        x0=x0,
        args=(
            y_measured,
            measured_spectrum_params,
            simulated_spectrum_params,
            x0_meaning,
        ),
        method="Nelder-Mead",
        bounds=bounds,
    )

    _update_params(
        measured_spectrum_params, simulated_spectrum_params, res.x, x0_meaning
    )

    if res.success:
        params.meas = measured_spectrum_params
        params.sim = simulated_spectrum_params

    return params, res


def _make_default_fit_settings(initial_guess: Params) -> FitSettings:
    """
    Create default fit settings based on the initial guess parameters.

    :param initial_guess: Initial parameters for the fit.
    :return: Default FitSettings object.
    """
    default_wav_shift_tolerance = 0.002

    meas = MeasuredSpectrumParamFitSettings(
        wav_start_min=initial_guess.meas.wav_start - default_wav_shift_tolerance,
        wav_start_max=initial_guess.meas.wav_start + default_wav_shift_tolerance,
        wav_start_fixed=False,
        wav_step_min=initial_guess.meas.wav_step * 0.9,
        wav_step_max=initial_guess.meas.wav_step * 1.1,
        wav_step_fixed=False,
        wav_2nd_min=-initial_guess.meas.wav_step * 1e03,
        wav_2nd_max=initial_guess.meas.wav_step * 1e03,
        wav_2nd_fixed=False,
    )

    sim = SimulatedSpectrumParamFitSettings(
        slitf_gauss_min=initial_guess.sim.slitf_gauss * 1e-3,
        slitf_gauss_max=initial_guess.sim.slitf_gauss * 1e3,
        slitf_gauss_fixed=False,
        slitf_lorentz_min=initial_guess.sim.slitf_lorentz * 1e-3,
        slitf_lorentz_max=initial_guess.sim.slitf_lorentz * 1e3,
        slitf_lorentz_fixed=False,
        baseline_min=-np.inf,
        baseline_max=np.inf,
        baseline_fixed=False,
        baseline_slope_min=0,
        baseline_slope_max=0,
        baseline_slope_fixed=True,
    )

    species = {
        species: SpectralSystemParamFitSettings(
            Tvib_min=250,
            Tvib_max=5e4,
            Tvib_fixed=False,
            Trot_min=250,
            Trot_max=5e4,
            Trot_fixed=False,
            intensity_min=0,
            intensity_max=np.inf,
            intensity_fixed=False,
        )
        for species in initial_guess.sim.species_params.keys()
    }

    settings = FitSettings(meas=meas, sim=sim, species=species)

    return settings


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
    """
    modifies the measured and simulated spectrum parameters in place
    """
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

    # x0 = [3000, 3000, 1]
    # x0_meaning = ["OHAX|Tvib", "OHAX|Trot", "OHAX|intensity"]

    # res = minimize(
    #     fun=_min_func,
    #     x0=x0,
    #     args=(
    #         y_measured,
    #         initial_measured_params,
    #         initial_simulated_params,
    #         x0_meaning,
    #     ),
    #     method="Nelder-Mead",
    #     # options={'disp': True}
    # )

    initial_guess = Params(
        meas=initial_measured_params,
        sim=initial_simulated_params,
    )

    params, res = fit(
        y_measured=y_measured,
        initial_guess=initial_guess,
        settings=None,
    )

    print(res)
