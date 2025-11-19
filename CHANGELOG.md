# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # "## [unreleased] - yyyy-mm-dd"


## [0.4.2] - 2025-11-18

### Added
- Added support for `MissingRCFBehaviour`: `NEAREST`, to return the closest wavelength's RCF of the ones that
  are missing. It's done so the actual photometer channel center can be used for the rest of calculations.
- Updated RCF coefficients with higher resolution ones.

### Changed
- Set default value of `per_nm` to True in the multiple places it appears.

### Fixed
- Correctly multiplying reflectance matrix by RCF, fixing wrong dimension match-up.

### Removed
- Removed support for `MissingRCFBehaviour`: `INTERPOLATE`, as interpolating the RCF does not make sense
  from a scientific point of view.


## [0.4.1] - 2025-11-17

### Added

- Implemented `ESICalculatorCustom` which allowes the user to load any ESI spectrum they want to use,
  and linearly interpolate for the desired wavelengths.

### Changed
- Set default value of `apply_correction` in `elref.get_reflectance` to `False` (previously `True`),
  to be aligned with the default values in `ELISettings`.

### Fixed
- Removed `slots=True` from `ELISettings` dataclass annotation to restore compatibility with Python 3.8 and 3.9.


## [0.4.0] - 2025-11-16

### Added
- Added support for configurable handling of missing RCF wavelengths by introducing the `MissingRCFBehavior`
  enum (ERROR, WARN, IGNORE, INTERPOLATE) as a new `missing_rcf` parameter in `ELISettings`.

### Changed
- Irradiance interpolation for wavelengths without RCFs is now performed by linearly
  interpolating the valid RCF values, instead of their coefficients.
- Unified all public ELI retrieval functions into a single `get_irradiance` entry point,
  replacing the former `get_eli`, `get_eli_bypass`, and `get_eli_from_extra_kernels`.
  The new function automatically handles the three geometry input modes
  (`MoonDatas`, `EarthPoint` + kernels, or extra observer kernels).
- Renamed and extended `get_interpolated_reflectance` to `get_reflectance`,
  providing the same unified geometry adaptability as `get_irradiance`.


### Fixed
- RCF is calculated using the moon phase angle instead of using the absolute moon phase angle, which was wrong.

### Removed

- Removed the option to run RIMO with interpolated coefficients, as this approach was incorrect.
  Consequently, the `interpolate_rolo_coefficients` parameter was removed from `ELISettings`.

## [0.3.0] - 2025-11-06

### Added
- Complete vectorisation based on wavelengths and lunar geometries.
- Class-level, thread-safe file cache in `ESICalculatorWehrli` for Wehrli data,
  improving performance and eliminating global state.
- Added `LGPLv3` license.

### Changed
- Removed default values that created objects in function definitions. Changed to None, with instances now created inside the function.
- Merged former `per_nm`-specific functions (`get_eli_bypass_per_nm` , `get_eli_per_nm`, `get_eli_per_nm_from_extra_kernels`)
  into a more unified `get_eli`interface.
  - The desired behaviour is now controlled via the `per_nm` attribute in `ELISettings`.
- `get_eli` functions now only accept iterables of wavelengths as input (single-wavelength calls removed).
- `ELICalculator` converted into an abstract base class, implemented by `ELICalculatorWehrli`.


## [0.2.1] - 2025-11-05

Initial version that serves as the baseline for tracking changes in the change log.


[unreleased]: https://github.com/GOA-UVa/rimopy/compare/v0.4.2...HEAD
[0.4.2]: https://github.com/GOA-UVa/rimopy/compare/v0.4.1...0.4.2
[0.4.1]: https://github.com/GOA-UVa/rimopy/compare/v0.4.0...0.4.1
[0.4.0]: https://github.com/GOA-UVa/rimopy/compare/v0.3.0...0.4.0
[0.3.0]: https://github.com/GOA-UVa/rimopy/compare/v0.2.1...0.3.0
[0.2.1]: https://github.com/GOA-UVa/rimopy/releases/tag/v0.2.1
