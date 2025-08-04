from typing import Dict

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from .enums import SpectralSystemEnum


class SpectralSystemParams(BaseModel):
    Tvib: float
    Trot: float
    intensity: float


class SpectralSystemParamFitSettings(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    Tvib_min: float = Field(default=0, validation_alias="TvibMin")
    Tvib_max: float = Field(default=np.inf, validation_alias="TvibMax")
    Tvib_fixed: bool = Field(validation_alias="TvibFixed")

    Trot_min: float = Field(default=0, validation_alias="TrotMin")
    Trot_max: float = Field(default=np.inf, validation_alias="TrotMax")
    Trot_fixed: bool = Field(validation_alias="TrotFixed")

    intensity_min: float = Field(default=0, validation_alias="intensityMin")
    intensity_max: float = Field(default=np.inf, validation_alias="intensityMax")
    intensity_fixed: bool = Field(validation_alias="intensityFixed")


class SimulatedSpectrumParams(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    slitf_gauss: float = Field(validation_alias="slitfGauss")
    slitf_lorentz: float = Field(validation_alias="slitfLorentz")
    baseline: float
    baseline_slope: float = Field(validation_alias="baselineSlope")
    species_params: Dict[SpectralSystemEnum, SpectralSystemParams] = Field(
        validation_alias="speciesParams"
    )


class SimulatedSpectrumParamFitSettings(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    slitf_gauss_min: float = Field(default=0, validation_alias="slitfGaussMin")
    slitf_gauss_max: float = Field(default=np.inf, validation_alias="slitfGaussMax")
    slitf_gauss_fixed: bool = Field(validation_alias="slitfGaussFixed")

    slitf_lorentz_min: float = Field(default=0, validation_alias="slitfLorentzMin")
    slitf_lorentz_max: float = Field(default=np.inf, validation_alias="slitfLorentzMax")
    slitf_lorentz_fixed: bool = Field(validation_alias="slitfLorentzFixed")

    baseline_min: float = Field(default=-np.inf, validation_alias="baselineMin")
    baseline_max: float = Field(default=np.inf, validation_alias="baselineMax")
    baseline_fixed: bool = Field(validation_alias="baselineFixed")

    baseline_slope_min: float = Field(
        default=-np.inf, validation_alias="baselineSlopeMin"
    )
    baseline_slope_max: float = Field(
        default=np.inf, validation_alias="baselineSlopeMax"
    )
    baseline_slope_fixed: bool = Field(validation_alias="baselineSlopeFixed")


class XAbsoluteParams(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    wav_start: float = Field(validation_alias="wavStart")
    wav_step: float = Field(validation_alias="wavStep")
    wav_2nd: float = Field(validation_alias="wav2nd")


class XShiftedParams(BaseModel):
    x: list[float]
    shift: float


class XAbsoluteParamFitSettings(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    wav_start_min: float = Field(default=0, validation_alias="wavStartMin")
    wav_start_max: float = Field(default=np.inf, validation_alias="wavStartMax")
    wav_start_fixed: bool = Field(validation_alias="wavStartFixed")

    wav_step_min: float = Field(default=0, validation_alias="wavStepMin")
    wav_step_max: float = Field(default=np.inf, validation_alias="wavStepMax")
    wav_step_fixed: bool = Field(validation_alias="wavStepFixed")

    wav_2nd_min: float = Field(default=-np.inf, validation_alias="wav2ndMin")
    wav_2nd_max: float = Field(default=np.inf, validation_alias="wav2ndMax")
    wav_2nd_fixed: bool = Field(validation_alias="wav2ndFixed")


class XShiftedParamFitSettings(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    shift_min: float = Field(default=-np.inf, validation_alias="shiftMin")
    shift_max: float = Field(default=np.inf, validation_alias="shiftMax")
    shift_fixed: bool = Field(validation_alias="shiftFixed")


class Params(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    x_params: XAbsoluteParams | XShiftedParams = Field(validation_alias="xParams")
    sim: SimulatedSpectrumParams


class FitSettings(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    x_params: XAbsoluteParamFitSettings | XShiftedParamFitSettings = Field(
        validation_alias="xParams"
    )
    sim: SimulatedSpectrumParamFitSettings
    species: Dict[SpectralSystemEnum, SpectralSystemParamFitSettings]
