# pyrospec

A Python package for simulating spectra of diatomic molecules and for molecular pyrometry of plasma.

It is based on the [massiveOES](https://bitbucket.org/OES_muni/massiveoes/src/master/README.md) project, but this package has no user interface. It contains the code needed for calculations and spectral databases.

## Installation

```bash
pip install pyrospec
```

## Quick Start

```python
from pyrospec import (
    get_spectrum, 
    generate_spectrum, 
    estimate_intensity, 
    fit,
    SpectralSystemEnum,
    Params,
    XAbsoluteParams,
    SimulatedSpectrumParams,
    SpectralSystemParams
)
```

## Main Functions

### `get_spectrum` - Single molecular system simulation

Simulate spectral line heights for a single molecular system at given temperatures.

```python
from pyrospec import get_spectrum, SpectralSystemEnum

# Simulate OH (A-X) spectrum
spectrum = get_spectrum(
    system=SpectralSystemEnum.OHAX,
    Trot=3000,  # Rotational temperature in K
    Tvib=4000,  # Vibrational temperature in K
    wmin=306,   # Minimum wavelength in nm
    wmax=320,   # Maximum wavelength in nm
    step=0.02   # Wavelength step in nm
)
```

### `generate_spectrum` - Multi-system spectrum generation

Generate a complete spectrum with multiple overlapping molecular systems, each with its own temperatures.

```python
from pyrospec import generate_spectrum, SimulatedSpectrumParams, SpectralSystemParams

# Define parameters for multiple species
species_params = {
    SpectralSystemEnum.OHAX: SpectralSystemParams(
        Tvib=4000, Trot=3000, intensity=1.0
    ),
    SpectralSystemEnum.C2SWAN: SpectralSystemParams(
        Tvib=3500, Trot=2800, intensity=0.5
    )
}

sim_params = SimulatedSpectrumParams(
    slitf_gauss=0.1,      # Gaussian slit function width
    slitf_lorentz=0.05,   # Lorentzian slit function width
    baseline=0.1,         # Baseline offset
    baseline_slope=0.0,   # Baseline slope
    species_params=species_params
)

spectrum = generate_spectrum(
    simspec_params=sim_params,
    wmin=300,
    wmax=330,
    step=0.02
)
```

### `estimate_intensity` - Initial intensity estimation

Get a reasonable initial guess for intensity when fitting experimental data.

```python
from pyrospec import estimate_intensity

# Estimate intensity for OH system
estimated_intensity = estimate_intensity(
    initial_guess=initial_params,  # Your initial parameter guess
    species=SpectralSystemEnum.OHAX,
    y_measured=measured_spectrum_data  # Your experimental data
)
```

### `fit` - Temperature fitting

Determine rotational and vibrational temperatures that best fit measured data.

```python
from pyrospec import fit, XAbsoluteParams

# Define initial guess
x_params = XAbsoluteParams(wav_start=306.0, wav_step=0.02, wav_2nd=0.0)
initial_guess = Params(x_params=x_params, sim=sim_params)

# Perform the fit
fitted_params, optimization_result = fit(
    y_measured=experimental_data,
    initial_guess=initial_guess,
    settings=None  # Use default fit settings
)

print(f"Fit successful: {optimization_result.success}")
print(f"Fitted temperatures: {fitted_params}")
```

## Available spectral databases

 * C₂ (d-a) aka Swan system
 * N₂ (C-B) aka second positive system
 * N₂⁺ (B-X) aka first negative system
 * NH (A-X)
 * NO (B-X)
 * OH (A-X)

None of the databases has a complete set of spectral lines, but it covers the most popular spectral regions used for molecular pyrometry in plasma.

## Scientific paper

* Voráč, J., Synek, P., Potočňáková, L., Hnilica, J., & Kudrle, V. (2017). Batch processing of overlapping molecular spectra as a tool for spatio-temporal diagnostics of power modulated microwave plasma jet. [Plasma Sources Science and Technology, 26(2), 025010.](https://doi.org/10.1088/1361-6595/aa51f0)
