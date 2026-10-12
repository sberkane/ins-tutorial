# Tutorial on Aided Inertial Navigation Systems

MATLAB examples accompanying **Tutorial on Aided Inertial Navigation Systems:
A Modern Treatment Using Lie-Group Theoretical Methods**, by Soulaimane Berkane.

The examples cover:

- Exact and classical IMU discretization.
- Invariant log-linear error propagation.
- Navigation error distributions.
- Uncertainty propagation on SE₂(3).
- GPS-aided navigation with the multiplicative and invariant EKFs, with and
  without magnetometer aiding.

## Requirements

MATLAB. The code was developed with R2026a and requires no additional toolboxes
or external datasets.

## Getting started

Download or clone this repository and open its folder in MATLAB. Run:

```matlab
reproduce_all
```

This runs all examples and saves the figures and numerical results in `results/`.
The example functions are in `Examples/`, and shared functions are in `Tools/`.

## Citation

Soulaimane Berkane, *Tutorial on Aided Inertial Navigation Systems:
Reproduction Code*, [version 1.0.0](https://github.com/sberkane/ins-tutorial/releases/tag/v1.0.0), 2026.

Citation metadata is available in [CITATION.cff](CITATION.cff).

## License

Copyright (c) 2026 Soulaimane Berkane. Released under the [MIT license](LICENSE).
