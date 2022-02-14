# Rimopy

![Version 0.0.26](https://img.shields.io/badge/version-0.0.26-informational)

Rimopy is a package consisting of an implementation of the ROLO model, following RIMO's
implementation, made in python.

The ROLO model comes from the RObotic Lunar Observatory, Kieffer and Stone, 2005.

The RIMO implementation is from multiple papers:
- Barreto et al., 2019: Evaluation of night-time aerosols measurements and lunar irradiance
models in the frame of the first multi-instrument nocturnal intercomparison campaign.
- Roman et al., 2020: Correction of a lunar-irradiance model for aerosol optical depth
retrieval and comparison with a star photometer.

## Requirements

- numpy >= 1.22.0
- scipy >= 1.7.3
- spiceypy >= 4.0.3

## Installation

```sh
pip install rimopy
```

### Kernels

In order to use the package, a directory with all the kernels must be downloaded.

That directory must contain the same elements as **/tests/kernels**, and the execution must
be allowed to read and write from that directory.

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
without applying the RIMO correction factor, and with a different method for the solar irradiance interpolation,
something like the following code block would work:

```python
from rimopy import esi

wavelength = 500
calc = esi.ESICalculator(esi.WehrliFile.ORIGINAL_WEHRLI, esi.ESIMethod.LINEAR_INTERPOLATION)
eli_settings = eli.ELISettings(False, False, True)
result = eli.get_eli_per_nm(wavelength, full_moon_Valladolid, kernels_path, calc, eli_settings)
```

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
python3 -m venv .venv
source .venv/bin/activate
```

The packages mentioned in [Requirements](#requirements) should be installed,
and finally the latest build should be installed.

```sh
pip install ../dist/rimopy-<latest_version>-py3-none-any.whl
```

### Executing the tests

At the moment only one test file exists, test_eli_aemet.py, but that will be different in the future.