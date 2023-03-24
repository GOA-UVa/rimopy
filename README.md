# Rimopy

![Version 0.1.8](https://img.shields.io/badge/version-0.1.8-informational)

Rimopy is a package consisting of an implementation of the ROLO model, following RIMO's
implementation, made in python.

The ROLO model comes from the RObotic Lunar Observatory, from the paper:
- Kieffer and Stone, 2005: The spectral irradiance of the Moon.

The RIMO implementation is from multiple papers:
- Barreto et al., 2019: Evaluation of night-time aerosols measurements and lunar irradiance
models in the frame of the first multi-instrument nocturnal intercomparison campaign.
- Roman et al., 2020: Correction of a lunar-irradiance model for aerosol optical depth
retrieval and comparison with a star photometer.

## Requirements

- numpy >= 1.22.0
- spiceypy >= 5.0.0
- spicedmoon >= 1.0.3

## Installation

```sh
pip install rimopy
```

### Kernels

In order to use the package, a directory with all the kernels must be downloaded.

That directory must contain the same elements as **/tests/kernels**, and the execution must
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

## Usage

The main functions are eli.get_eli and eli.get_eli_per_nm. Both return the extraterrestrial
lunar irradiance for a set of wavelengths or for only one, in nanometers. This irradiance is
calculated at a concrete location on Earth's surface, and time instance.
The first one returns the value in Wm⁻², and the second one in Wm⁻²/nm.

In order to calculate the extraterrestrial lunar irradiance per nm for 500 nm at the city of
Valladolid, the morning of the 2022-01-17, one would do something like the following code block:

```python
from rimopy import eli

wavelength = 500
full_moon_Valladolid = eli.EarthPoint(41.652251, -4.7245321, "2022-01-17 03:00:00", 700)
kernels_path = "./kernels"
result = eli.get_eli_per_nm(wavelength, full_moon_Valladolid, kernels_path)
```

In order to calculate the extraterrestrial lunar irradiance for a set of wavelengths at the
same time and space conditions, the next code block would work:

```python
wavelengths = [500, 510, 520, 530, 540, 550]
results = eli.get_eli(wavelengths, full_moon_Valladolid, kernels_path)
```

These calculations can be customized, defining the settings and the methods used for the calculation of
the extraterrestrial solar irradiance. For example, if someone wanted to calculate the lunar irradiance
applying the RIMO correction factor, and with a different source data for the extraterrestrial solar irradiance,
in this case the original Wehrli 1985 data, something like the following code block would work:

```python
from rimopy import esi

wavelength = 500
calc = esi.ESICalculator(esi.WehrliFile.ORIGINAL_WEHRLI, esi.ESIMethod.LINEAR_INTERPOLATION)
eli_settings = eli.ELISettings(True, False, True)
result = eli.get_eli_per_nm(wavelength, full_moon_Valladolid, kernels_path, calc, eli_settings)
```

### Main functions

The most important functions are the ones that obtain the extraterrestrial lunar irradiance, so they
are all contained in the **eli** module.

* **get_eli_bypass** - returns the expected extraterrestrial lunar irradiation of a wavelength for
any observer/solar selenographic coordinates, in Wm⁻².
* **get_eli** - returns the expected extraterrestrial lunar irradiation of a wavelength in any
geographic coordinates, in Wm⁻².
* **get_eli_from_extra_kernels** - returns the expected extraterrestrial lunar irradiation of
a wavelength in any geographic coordinates, in Wm⁻², using data from extra kernels for the
observer body.
* **get_eli_bypass_per_nm** - returns the expected extraterrestrial lunar irradiation of a
wavelength for any observer/solar selenographic coordinates, in Wm⁻²/nm.
* **get_eli_per_nm** - returns the expected extraterrestrial lunar irradiation of a wavelength in
any geographic coordinates, in Wm⁻²/nm.
* **get_eli_per_nm_from_extra_kernels** - returns the expected extraterrestrial lunar irradiation of
a wavelength in any geographic coordinates, in Wm⁻²/nm, using data from extra kernels for the
observer body.

### Configuration options

**ELISettings**

ELISettings is a dataclass that contains settings that will modify the methodology of calculating the ELI.

These settings are:
- **apply_correction**: If True the result will have been multiplied by the RCF (Rimo Correction Factor),
which corrects the data for the photometers' calibration. Otherwise it won't. The default value is *False*.
- **interpolate_rolo_coefficients**: If True the reflectance will be calculated linearly interpolating
the ROLO coefficients. Otherwise it will be calculated interpolating the surrounding reflectances,
calculated with empirical coefficients. The default value is *False*.
- **adjust_apollo**: If True the ROLO model reflectance will be adjusted using Apollo spectra, in case
it's calculated interpolating surrounding reflectances. The Apollo spectra is the spectra generated
with the Moon samples from Apollo 16th mission. The default value is *True*.

**ESICalculator**

ESICalculator is the class of which instances contain the functions that calculate the extraterrestrial
solar irradiance. This calculations will vary depending on the values of ESICalculator attributes.

These attributes are:
- **wehrli_file**: Wehrli data source that will be used in the calculation of the ESI. It could be the original
data or some filtered data. The default value is "ASC_WEHRLI", the Wehrli data passed throught some filters, and it's
the one used in AEMET's RimoApp and others. Although it's the default value, it might not be the best for obtaining
"Wm⁻²" data (**get_eli()**), as in this file the solar irradiance in "Wm⁻²" is obtained from the "Wm⁻²/nm" data, instead of
having passed the original "Wm⁻²" data through the same filters.
- **method**: Interpolation method that will be used in the calculation of the ESI. The default one is "LINEAR_INTERPOLATION".
- **gaussian_filter_params**: Parameters for the gaussian filter method, in case that that one is the chosen one. It contains
the radius and the standard deviation, which by default both are equal to one.

## Structure

The package is divided in multiple modules, each dealing with different calculations and
functionalities:

- **eli**: Main module, which calculates the Extraterrestrial Lunar Irradiance for a given data.
- **esi**: Calculates the Extraterrestrial Solar Irradiance using data from Wehrli 1985. The methodology for calculating this data
can be chosen by the user, creating an instance of ESICalculator selecting the methodology and data source they want.
- **coefficients**: Contains the coefficient data from the ROLO model.
- **correction_factor**: Calculates the correction factor as stated in RIMO papers.
- **spice_iface**: Encapsulates the access to functionalities from SPICE Toolbox.
- **types**: Contains types, like MoonData class, which represents some of the needed data
for the calculation of extraterrestrial lunar irradiance.

![ModuleStructure UML Diagram](docs/ModuleStructure.png)

## Build

```sh
python3 -m build
```

## Testing

### Testing setup

A virtual environment should be created inside the 'tests' directory.

```sh
cd tests
python3 -m venv .venv
source .venv/bin/activate
```

The packages mentioned in [Requirements](#requirements) should be installed,
and finally the latest build should be installed.

```sh
pip install ../dist/rimopy-<latest_version>-py3-none-any.whl
```

### Executing the tests

The tests shall be executed from the tests directory, either directly or using pytest.

```sh
# Directly
./test_eli_aemet.py

# Using pytest
# pip install pytest
pytest -v
```
At the moment there is only one test module, test_eli_aemet.py, which tests rimopy against output
from AEMET's RimoApp.