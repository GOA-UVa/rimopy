# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # "## [unreleased] - yyyy-mm-dd"

## [unreleased] - yyyy-mm-dd

### Changed
- Removed default values that created objects in function definitions. Changed to None and the instance is created inside the function.
- Removed `per_nm` specific functions, and added it as an attribute of `ELISettings`
- `get_eli` functions now only accept an iterable of wavelengths as input, removing the option for a sole wavelength.

## [0.2.1] - 2025-11-05

Initial version that serves as the baseline for tracking changes in the change log.

[unreleased]: https://gitlab.com/GOA-UVa/rimopy/compare/v0.2.1...HEAD
[0.2.1]: https://gitlab.com/GOA-UVa/rimopy/-/tags/v0.2.1
