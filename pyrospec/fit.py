import copy

import numpy as np
from .enums import SpectralSystemEnum
from .models import (
    FitSettings,
    XAbsoluteParamFitSettings,
    XShiftedParamFitSettings,
    XAbsoluteParams,
    XShiftedParams,
    Params,
    SimulatedSpectrumParamFitSettings,
    SimulatedSpectrumParams,
    SpectralSystemParamFitSettings,
)
from scipy.optimize import minimize
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

    x_params = params.x_params
    simulated_spectrum_params = params.sim

    x0 = []
    x0_meaning = []

    if settings is None:
        settings = _make_default_fit_settings(initial_guess)

    bounds = []

    for param_name, param_value in x_params.__dict__.items():
        if (
            hasattr(settings.x_params, f"{param_name}_fixed")
            and getattr(settings.x_params, f"{param_name}_fixed")
            or param_name == "x"
        ):
            # Only not fixed parameters go into the optimization
            continue
        x0.append(param_value)
        x0_meaning.append(param_name)

        lb = getattr(settings.x_params, f"{param_name}_min", -np.inf)
        ub = getattr(settings.x_params, f"{param_name}_max", np.inf)
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
            x_params,
            simulated_spectrum_params,
            x0_meaning,
        ),
        method="Nelder-Mead",
        bounds=bounds,
        options={"maxiter": 10000},
    )

    _update_params(x_params, simulated_spectrum_params, res.x, x0_meaning)

    if res.success:
        params.x_params = x_params
        params.sim = simulated_spectrum_params

    return params, res


def _make_default_fit_settings(initial_guess: Params) -> FitSettings:
    """
    Create default fit settings based on the initial guess parameters.

    :param initial_guess: Initial parameters for the fit.
    :return: Default FitSettings object.
    """
    default_wav_shift_tolerance = 0.002

    if isinstance(initial_guess.x_params, XAbsoluteParams):
        x_params = XAbsoluteParamFitSettings(
            wav_start_min=initial_guess.x_params.wav_start
            - default_wav_shift_tolerance,
            wav_start_max=initial_guess.x_params.wav_start
            + default_wav_shift_tolerance,
            wav_start_fixed=False,
            wav_step_min=initial_guess.x_params.wav_step * 0.9,
            wav_step_max=initial_guess.x_params.wav_step * 1.1,
            wav_step_fixed=False,
            wav_2nd_min=-initial_guess.x_params.wav_step * 1e03,
            wav_2nd_max=initial_guess.x_params.wav_step * 1e03,
            wav_2nd_fixed=False,
        )
    elif isinstance(initial_guess.x_params, XShiftedParams):
        x_params = XShiftedParamFitSettings(
            shift_min=-default_wav_shift_tolerance,
            shift_max=default_wav_shift_tolerance,
            shift_fixed=False,
        )
    else:
        raise ValueError("Unsupported x_params type in initial guess.")

    slitf_sum = initial_guess.sim.slitf_gauss + initial_guess.sim.slitf_lorentz

    sim = SimulatedSpectrumParamFitSettings(
        slitf_gauss_min=slitf_sum * 1e-3,
        slitf_gauss_max=slitf_sum * 3,
        slitf_gauss_fixed=False,
        slitf_lorentz_min=slitf_sum * 1e-3,
        slitf_lorentz_max=slitf_sum * 3,
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

    settings = FitSettings(x_params=x_params, sim=sim, species=species)

    return settings


def get_residuals(
    measured: np.ndarray,
    x_params: XAbsoluteParams | XShiftedParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    simulations: dict[SpectralSystemEnum, SpecDB],
) -> np.ndarray:

    x_measured, step, wmin = get_x_measured(x_params, len(measured))

    simulated_spec: spectrum.Spectrum = generate_spectrum(
        simspec_params=simulated_spectrum_params,
        step=step,
        wmin=wmin,
        wmax=np.max(x_measured),
        simulations=simulations,
    )

    measured_spec = spectrum.Spectrum(x=x_measured, y=measured)

    sim_spec_matched, exp_spec = spectrum.match_spectra(
        sim_spec=simulated_spec, exp_spec=measured_spec
    )

    return sim_spec_matched.y - exp_spec.y


def get_x_measured(
    x_params: XAbsoluteParams | XShiftedParams, length: int
) -> tuple[np.ndarray, float, float]:
    if isinstance(x_params, XAbsoluteParams):
        x_axis = (
            x_params.wav_start
            + np.arange(length, dtype=float) * x_params.wav_step
            + x_params.wav_2nd * np.arange(length, dtype=float) ** 2
        )
        step = abs(x_params.wav_step)

    elif isinstance(x_params, XShiftedParams):
        x_axis = np.array(x_params.x, dtype=float) + x_params.shift
        step = float(np.mean(np.diff(x_axis)))

    wmin = float(np.min(x_axis))
    return x_axis, step, wmin


def _min_func(
    x0,
    y_measured,
    x_params: XAbsoluteParams | XShiftedParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    x0_meaning,
):
    _update_params(x_params, simulated_spectrum_params, x0, x0_meaning)
    residuals = get_residuals(
        measured=y_measured,
        x_params=x_params,
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
    x_params: XAbsoluteParams | XShiftedParams,
    simulated_spectrum_params: SimulatedSpectrumParams,
    x0,
    x0_meaning,
):
    """
    modifies the measured and simulated spectrum parameters in place
    """
    species_param_keywords = ["Tvib", "Trot", "intensity"]
    for param, value in zip(x0_meaning, x0):
        if param in x_params.__dict__:
            x_params.__dict__[param] = value
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


def estimate_intensity(
    initial_guess: Params, species: SpectralSystemEnum, y_measured: np.ndarray
) -> float:
    """
    Estimate a better initial guess for intensity based on the initial guess parameters.
    """
    x_axis, step, wmin = get_x_measured(initial_guess.x_params, len(y_measured))
    simulated_spec = generate_spectrum(
        simspec_params=initial_guess.sim,
        step=step,
        wmin=wmin,
        wmax=np.max(x_axis),
        simulations={species: SpecDB(species.filename)},
    )
    measured_spec = spectrum.Spectrum(x=x_axis, y=y_measured)
    sim_spec_matched, exp_spec = spectrum.match_spectra(
        sim_spec=simulated_spec, exp_spec=measured_spec
    )
    est_intensity = np.max(exp_spec.y) / np.max(sim_spec_matched.y)
    if np.isnan(est_intensity) or est_intensity < 0 or np.isinf(est_intensity):
        est_intensity = 1
    return est_intensity
