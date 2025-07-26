from typing import Dict

import numpy as np
from pydantic import BaseModel

from .enums import SpectralSystemEnum


class SpectralSystemParams(BaseModel):
    Tvib: float
    Trot: float
    intensity: float


class SpectralSystemParamFitSettings(BaseModel):
    Tvib_min: float = 0
    Tvib_max: float = np.inf
    Tvib_fixed: bool

    Trot_min: float = 0
    Trot_max: float = np.inf
    Trot_fixed: bool

    intensity_min: float = 0
    intensity_max: float = np.inf
    intensity_fixed: bool


class SimulatedSpectrumParams(BaseModel):
    slitf_gauss: float
    slitf_lorentz: float
    baseline: float
    baseline_slope: float
    species_params: Dict[SpectralSystemEnum, SpectralSystemParams]


class SimulatedSpectrumParamFitSettings(BaseModel):
    slitf_gauss_min: float = 0
    slitf_gauss_max: float = np.inf
    slitf_gauss_fixed: bool

    slitf_lorentz_min: float = 0
    slitf_lorentz_max: float = np.inf
    slitf_lorentz_fixed: bool

    baseline_min: float = -np.inf
    baseline_max: float = np.inf
    baseline_fixed: bool

    baseline_slope_min: float = -np.inf
    baseline_slope_max: float = np.inf
    baseline_slope_fixed: bool


class MeasuredSpectrumParams(BaseModel):
    wav_start: float
    wav_step: float
    wav_2nd: float


class MeasuredSpectrumParamFitSettings(BaseModel):
    wav_start_min: float = 0
    wav_start_max: float = np.inf
    wav_start_fixed: bool

    wav_step_min: float = 0
    wav_step_max: float = np.inf
    wav_step_fixed: bool

    wav_2nd_min: float = -np.inf
    wav_2nd_max: float = np.inf
    wav_2nd_fixed: bool


class Params(BaseModel):
    meas: MeasuredSpectrumParams
    sim: SimulatedSpectrumParams


class FitSettings(BaseModel):
    meas: MeasuredSpectrumParamFitSettings
    sim: SimulatedSpectrumParamFitSettings
    species: Dict[SpectralSystemEnum, SpectralSystemParamFitSettings]
