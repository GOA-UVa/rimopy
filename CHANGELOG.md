# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # "## [unreleased] - yyyy-mm-dd"

## [unreleased] - yyyy-mm-dd

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


[unreleased]: https://gitlab.com/GOA-UVa/rimopy/compare/v0.2.1...HEAD
[0.2.1]: https://gitlab.com/GOA-UVa/rimopy/-/tags/v0.2.1
