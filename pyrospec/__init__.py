# Version and contact info for pyrospec
__version__ = "0.1.0"
__author__ = "Jan Vorac"
__email__ = "janvorac@proton.me"

# Direct imports for user convenience
from .spectrum import Spectrum, match_spectra
from .fit import fit
from .enums import SpectralSystemEnum
from .specdata import SpecDB, generate_spectrum, list_available_databases
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
