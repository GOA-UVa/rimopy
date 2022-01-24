# Wehrli Gauss

This utility takes the data from 1985 Wehrli Standard Extraterrestrial Solar
Irradiance Spectrum (which can be found at
[nrel.gov/grid/solar-resource/spectra-wehrli.html](https://www.nrel.gov/grid/solar-resource/spectra-wehrli.html)
), and filters it.

It takes the original data, interpolates it so we have data for every 2 nanometers from 199.5 to 3000,
and filters it with a Gaussian 1D Filter.

The original data is in **./data/wehrli_original.csv** and the generated data is in
**./data/wehrli_filtered.csv**.

