[![Latest Tag][version-shield]][version-url]
[![License: LGPL 3.0][license-shield]][license-url]
[![Stargazers][stars-shield]][stars-url]

<br />
<div align="center">
  <a href="https://github.com/goa-uva/rimopy">
    <img src="https://raw.githubusercontent.com/GOA-UVa/rimopy/master/docs/images/rimo_logo.png" alt="RIMO logo" height="80" style="height: 80px !important;">
  </a>

  <h3 align="center">rimopy</h3>

  <p align="center">
    ROLO Implementation for Moon photometry Observation (RIMO), in Python.
  </p>
</div>

## About the project

`rimopy` is a Python package that implements the **RIMO** (ROLO Implementation for Moon photometry
Observation) lunar-irradiance model for night-time aerosol optical depth (AOD) retrieval from lunar
photometry (Barreto et Al., 2019). It also includes the RIMO Correction Factor (RCF), an empirical
adjustment that improves AOD retrieval accuracy (Román et al., 2020). RIMO builds on the ROLO
(RObotic Lunar Observatory) model of lunar spectral irradiance (Kieffer & Stone, 2005).

## Getting started

### Prerequisites

* python >= 3.8

### Kernels

In order to use the package, a directory with all the needed SPICE kernels must be downloaded.

That directory must contain the same elements as `rimopy/tests/kernels`, and the execution must
be allowed to read and write from that directory.

Alternatively, kernels can be downloaded manually from the following urls:
- [https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/de421.bsp](https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/de421.bsp)
- [https://naif.jpl.nasa.gov/pub/naif/pds/wgc/kernels/pck/earth_070425_370426_predict.bpc](https://naif.jpl.nasa.gov/pub/naif/pds/wgc/kernels/pck/earth_070425_370426_predict.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0011.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0011.tls)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc)

### Installation

You can install `rimopy` either from the source directory (for development) or directly from PyPI.

From source (editable mode):
```sh
pip install -e .
```

From PyPI (recommended):
```sh
pip install rimopy
```

## Usage

The main function is `eli.get_eli`. It returns the extraterrestrial
lunar irradiance for a set of wavelengths in nanometers. This irradiance is
calculated at a specific location on Earth's surface, for a set of date and time instants.
It returns the irradiance in Wm⁻² or in Wm⁻²/nm, depending on the configuration.

In order to calculate the extraterrestrial lunar irradiance per nm for 500 and 550 nm at the city of
Valladolid, the morning of the 2022-01-17, one would do something like the following code block:

```python
from rimopy import eli

wavelengths = [500, 550]
full_moon_Valladolid = eli.EarthPoint(41.652251, -4.7245321, "2022-01-17 03:00:00", 700)
kernels_path = "./kernels"
eli_settings = eli.ELISettings(per_nm=True)
result = eli.get_eli(wavelengths, full_moon_Valladolid, kernels_path, eli_settings=eli_settings)

# result:
# array([[3.05176826e-06],
#       [3.25669721e-06]])
```

These calculations can be customized, defining the settings and the methods used for the calculation of
the extraterrestrial solar irradiance. For example, if someone wanted to calculate the lunar irradiance per nm
applying the RIMO correction factor, and with a different source data for the extraterrestrial solar irradiance,
in this case the original Wehrli 1985 data, something like the following code block would do the work:
```python
from rimopy import esi

wavelengths = [500, 550]
calc = esi.ESICalculatorWehrli(esi.WehrliFile.ORIGINAL_WEHRLI, esi.ESIMethod.LINEAR_INTERPOLATION)
eli_settings = eli.ELISettings(True, False, True, True)
result = eli.get_eli(wavelengths, full_moon_Valladolid, kernels_path, calc, eli_settings)

# result:
# array([[3.29045832e-06],
#       [3.52460413e-06]])
```

### Main functions

The most important functions are the ones that obtain the extraterrestrial lunar irradiance, so they
are all contained in the **eli** module.

* **get_eli_bypass** - returns the expected extraterrestrial lunar irradiation of a wavelength for
any observer and solar selenographic coordinates.
* **get_eli** - returns the expected extraterrestrial lunar irradiation of a wavelength in any
geographic coordinates.
* **get_eli_from_extra_kernels** - returns the expected extraterrestrial lunar irradiation of
a wavelength in any geographic coordinates, using data from extra kernels for the
observer body.

### Configuration options

**ELISettings**

`ELISettings` is a dataclass that contains settings that will modify the methodology of calculating the ELI.

These settings are:
- **apply_correction**: If True the result will have been multiplied by the RCF (Rimo Correction Factor),
which corrects the data for the photometers' calibration. Otherwise it won't. The default value is *False*.
- **interpolate_rolo_coefficients**: If True the reflectance will be calculated linearly interpolating
the ROLO coefficients. Otherwise it will be calculated interpolating the surrounding reflectances,
calculated with empirical coefficients. The default value is *False*.
- **adjust_apollo**: If True the ROLO model reflectance will be adjusted using Apollo spectra, in case
it's calculated interpolating surrounding reflectances. The Apollo spectra is the spectra generated
with the Moon samples from Apollo 16th mission. The default value is *True*.
- **per_nm** : If True the ELI will be in Wm⁻²/nm, otherwise it will be in Wm⁻². Default is *False*.

**ESICalculator**

`ESICalculator` is the interface of which implementations contain the functions that calculate the extraterrestrial
solar irradiance. This calculations will vary depending on the implementation.

Currently there's only one implementation: `ESICalculatorWehrli`. Its attributes are:
- **wehrli_file**: Wehrli data source that will be used in the calculation of the ESI. It could be the original
data or some filtered data. The default value is "ASC_WEHRLI", the Wehrli data passed throught some filters, and it's
the one used in AEMET's RimoApp and others. Although it's the default value, it might not be the best for obtaining
"Wm⁻²" data (**get_eli()**), as in this file the solar irradiance in "Wm⁻²" is obtained from the "Wm⁻²/nm" data, instead of
having passed the original "Wm⁻²" data through the same filters.
- **method**: Interpolation method that will be used in the calculation of the ESI. The default one is "LINEAR_INTERPOLATION".
- **gaussian_filter_params**: Parameters for the gaussian filter method, in case that that one is the chosen one. It contains
the radius and the standard deviation, which by default both are equal to one.


## Roadmap

- [ ] Adding tests based on AEMET RimoApp output that don't require to compute the lunar geometries (have them precomputed
- [ ] Adding unit tests for every submodule.
- [ ] Adding a TSIS-based ESICalculator implementation.

## Structure

The package is divided in multiple modules, each dealing with different calculations and
functionalities:

- **eli**: Main module, which calculates the Extraterrestrial Lunar Irradiance for a given input.
- **elref**: Calculates the Extraterrestrial Lunar Reflectance following RIMO.
- **esi**: Calculates the Extraterrestrial Solar Irradiance using data from Wehrli 1985. The methodology for calculating this data
can be chosen by the user, creating an instance of ESICalculator selecting the methodology and data source they want.
- **coefficients**: Contains the coefficient data from the ROLO model.
- **correction_factor**: Calculates the correction factor as stated in RIMO papers.
- **spice_iface**: Encapsulates the access to functionalities from SPICE Toolbox.
- **types**: Contains types, like MoonData class, which represents some of the needed data
for the calculation of extraterrestrial lunar irradiance.

![ModuleStructure UML Diagram](https://raw.githubusercontent.com/GOA-UVa/rimopy/master/docs/ModuleStructure.png)

## Development Guide

### Setting up the environment

To set up the development environment, install the pre-commit hooks as follows:
```sh
pre-commit install
```
This ensures that code style checks and formatting rules are automatically applied before each commit.

### Build

Build the package with the following command:
```sh
python3 -m build
```
This command creates a distribution package under the dist/ directory, which can be uploaded to PyPI
or used for testing installations.

### Testing

Currently, the only tests available are the ones at `rimopy/tests/test_eli_aemet.py`, which
tests `rimopy` against output from AEMET's [RimoApp](https://testbed.aemet.es/rimoapp/).

RimoApp output test-suit can be run directly:
```sh
./rimopy/tests/test_eli_aemet.py
```

Or using `pytest`, if installed:
```sh
pytest -v
```

> Tests must be run in an environment where `rimopy` is installed or available.

## Authors

- [Javier Gatón Herguedas](mailto:gaton@goa.uva.es) - *Maintainer* - [GOA-UVa](https://goa.uva.es)

## License

Distributed under the LGPL-v3 License. See [LGPL v3](./LICENSE) for more information.

## References

- Barreto, Á., Román, R., Cuevas, E., Pérez-Ramírez, D., Berjón, A. J., Kouremeti, N., Kazadzis, S., Gröbner, J., Mazzola, M., Toledano, C., Benavent-Oltra, J. A., Doppler, L., Juryšek, J., Almansa, A. F., Victori, S., Maupin, F., Guirado-Fuentes, C., González, R., Vitale, V., Goloub, P., Blarel, L., Alados-Arboledas, L., Woolliams, E., Taylor, S., Antuña, J. C., & Yela, M. (2019). Evaluation of night-time aerosols measurements and lunar irradiance models in the frame of the first multi-instrument nocturnal intercomparison campaign. Atmospheric Environment, 202, 190–211. [https://doi.org/10.1016/j.atmosenv.2019.01.006](https://doi.org/10.1016/j.atmosenv.2019.01.006).
- Román, R., González, R., Toledano, C., Barreto, Á., Pérez-Ramírez, D., Benavent-Oltra, J. A., Olmo, F. J., Cachorro, V. E., Alados-Arboledas, L., & de Frutos, Á. M. (2020). Correction of a lunar-irradiance model for aerosol optical depth retrieval and comparison with a star photometer. *Atmospheric Measurement Techniques, 13*(11), 6293–6310. [https://doi.org/10.5194/amt-13-6293-2020](https://doi.org/10.5194/amt-13-6293-2020).
- Kieffer, H. H., & Stone, T. C. (2005). The spectral irradiance of the Moon. *The Astronomical Journal, 129*(6), 2887–2901. [https://doi.org/10.1086/430185](https://doi.org/10.1086/430185).


[stars-shield]: https://img.shields.io/github/stars/goa-uva/rimopy.svg?style=for-the-badge
[stars-url]: https://github.com/goa-uva/rimopy/stargazers
[version-shield]: https://img.shields.io/github/v/tag/goa-uva/rimopy?style=for-the-badge
[version-url]: https://github.com/goa-uva/rimopy/tags
[license-shield]: https://img.shields.io/github/license/goa-uva/rimopy.svg?style=for-the-badge
[license-url]: https://github.com/goa-uva/rimopy/blob/master/LICENSE
