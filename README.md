# Inertial navigation tutorial: numerical examples

MATLAB code for the numerical experiments in the inertial navigation tutorial.
The examples cover exact discretization, invariant error dynamics, uncertainty
representations, covariance propagation, and GPS-aided filtering on SE₂(3).

## Requirements

Validated with MATLAB R2026a on macOS (Apple silicon). No additional MathWorks toolboxes,
third-party packages, or downloaded datasets are required. GNU Octave has not
been validated. The examples generate vector PDF figures using MATLAB graphics.

## Reproduce the results

Open MATLAB in the repository root, then run:

```matlab
run_tests
reproduce_all
```

`run_tests` checks the numerical building blocks. `reproduce_all` runs all five
experiments, including the complete 50-seed filter comparison, and writes results
to `results/`. Existing files with the same names in that directory are replaced.

Each example can also be run independently:

```matlab
addpath('Examples')
result = loglinear_error;
result = compare_filters(fullfile(pwd, 'results'));
```

All example functions accept an optional output directory and return their
numerical results. They restore the caller's random-number generator, MATLAB
path, and graphics interpreter defaults on return.

## Figure map

Figure numbers refer to the October 11, 2026 tutorial PDF. Output filenames also match the figure references in the manuscript source.

| Figure | Example function | PDF output |
| --- | --- | --- |
| 6 | `discretization` | `exact_vs_classical.pdf` |
| 8 | `loglinear_error` | `loglinear_vs_mekf.pdf` |
| 9 | `error_distributions` | `distribution_growth.pdf` |
| 10 | `uncertainty_propagation` | `uncertainty_propagation.pdf` |
| 12 | `compare_filters` | `Example_circular.pdf` |

`compare_filters` also exports `Example_circular_log.pdf`, showing the full
60-second comparison on logarithmic error axes, and PNG copies of both figures.
The main `Example_circular.pdf` shows the first 10 seconds in its error panels;
its trajectory panels cover the full 60 seconds.

## Experiments and fixed settings

- **Discretization:** noise-free circular motion for 20 seconds, speed 5 m/s,
  radius 12.5 m, and sampling intervals of 0.01, 0.05, and 0.1 seconds.
- **Log-linear error:** straight and circular motion over 5 seconds at a
  0.01-second step, starting with a 20-degree attitude error about the normalized
  axis `[1; -1; 0.5]`. Both nonlinear trajectories receive the same IMU input.
- **Error distributions:** 20,000 tangent-space samples using seed 7. A fixed
  subset is displayed at nominal x positions from −30 to 30 m in 10 m steps.
- **Uncertainty propagation:** 500 Monte Carlo trajectories over 15 seconds at
  a 0.05-second step, using seed 7 and yaw-rate noise with standard deviation
  20 degrees/s per discrete sample. A small diagonal regularization permits
  sampling the singular process-noise covariance. The orange ellipse uses the
  sample covariance and is centered on the nominal endpoint.
- **Filter comparison:** a 60-second circle with radius 20 m and speed 3 m/s,
  IMU at 100 Hz, GPS at 10 Hz, and optional magnetometer aiding at 20 Hz.

Random experiments explicitly use MATLAB's `twister` generator. The filter
comparison uses every seed from 101 through 150. Its initial condition is a
30-degree pitch error, with exact initial velocity and position. This scenario
was selected in exploratory runs using seeds 1 through 5; the displayed results
use the separate 50-seed set.

For each seed, both filters receive the same measurements. The IMU and GPS noise
realizations are also paired between the configurations with and without
magnetometer aiding. The initial covariances represent the same local physical
uncertainty to first order. The MEKF uses classical mean propagation; the InvEKF
uses exact zero-order-held propagation. Neither filter applies covariance resets.

Error curves in Figure 12 are pointwise medians. Shaded bands are the 10th–90th
percentiles across realizations, not confidence intervals or predicted
covariance bounds. Trajectories use seed 101. The attitude error is the norm of
the rotation logarithm; velocity and position errors use Euclidean norms.

## Numerical outputs

Each of the first four examples saves a MAT file named after its function,
containing a `result` structure. The filter comparison saves:

- `paired_run_metrics.csv`: one row per seed, sensing configuration, and error
  quantity, with RMSE over 0–10 s, 0–60 s, and 50–60 s, and the peak error.
- `validated_results.mat`: parameters, ground truth, seeds, per-run metrics,
  error norms for all realizations, and the first realization of each scenario.

Attitude metrics in the CSV and MAT files are in radians; velocity metrics are
in m/s and position metrics are in metres. CSV quantity codes are 1 for attitude,
2 for velocity, and 3 for position; the magnetometer flag is 0 or 1.

## Repository layout

`Examples/` contains the five experiments and their private setup helper.
`Tools/` contains the shared SO(3), SE₂(3), motion, and propagation functions.
Generated figures and numerical data are written to `results/`.

These simulations illustrate the stated models and scenarios. They do not test
sensor biases, real-data performance, universal filter superiority, or
statistical covariance consistency. The log-linear example measures error-model
mismatch and does not perform measurement-based state estimation.

## Citation

Soulaimane Berkane, *Tutorial on Aided Inertial Navigation Systems:
Reproduction Code*, version 1.0.0, 2026.
https://github.com/sberkane/ins-tutorial/tree/v1.0.0

Citation metadata is provided in `CITATION.cff`. The accompanying tutorial is
*Tutorial on Aided Inertial Navigation Systems: A Modern Treatment Using
Lie-Group Theoretical Methods*, LaRSA Technical Report LaRSA-TR-2026-01.

## License

Copyright (c) 2026 Soulaimane Berkane. Released under the MIT license; see `LICENSE`.
