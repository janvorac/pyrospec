# pyrospec

A Python package for simulating spectra of diatomic molecules and for molecular pyrometry of plasma.

It is based on the [massiveOES](https://bitbucket.org/OES_muni/massiveoes/src/master/README.md) project, but this package has no user interface. It contains the code needed for calculations and spectral databases.

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
