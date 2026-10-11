# Rounded reproduction checks

These values are rounded exactly as in the manuscript and are not substitute
data. Small Monte Carlo differences are expected. A reproduced fit should agree
at the displayed precision and should pass the diagnostic thresholds in
`R/00_setup.R`.

## Primary effort psychophysics

Finite-sample participant-average MDD-minus-control contrasts:

| Quantity | Mean | 95% credible interval |
|---|---:|---:|
| Rating intercept | 0.187 | [0.149, 0.226] |
| Effort discriminability | -0.072 | [-0.110, -0.033] |

The public output also reports the hierarchical population contrast separately.

## Actual-execution sensitivities

| Adjustment | Discriminability contrast | 95% credible interval | pd |
|---|---:|---:|---:|
| Force AUC | -0.079 | [-0.150, -0.006] | .983 |
| Peak force | -0.126 | [-0.220, -0.030] | .995 |
| Mean force | -0.052 | [-0.135, 0.031] | .890 |

## Dynamic observer

| Quantity | MDD-minus-control | 95% credible interval | pd |
|---|---:|---:|---:|
| Prior uncertainty | 4.80 | [-0.48, 10.32] | .963 |
| Sensory uncertainty | 4.59 | [3.11, 6.17] | >.999 |
| Effective prior weighting | 0.030 | [0.004, 0.057] | .990 |

## Early force estimate

| Quantity | MDD-minus-control | 95% credible interval | pd |
|---|---:|---:|---:|
| Required-effort coefficient | -0.140 | [-0.232, -0.047] | .999 |
| Early-estimate coupling | 0.105 | [-0.030, 0.242] | .938 |
| Early-estimate by effort | -0.131 | [-0.195, -0.067] | >.999 |

Simple MDD-minus-control early-estimate coupling contrasts at 10, 30, 50, and
70% MVC are approximately 0.280, 0.161, 0.042, and -0.077, respectively.

## Decision making

| Quantity | MDD-minus-control | 95% credible interval |
|---|---:|---:|
| Effort sensitivity | -3.037 | [-4.637, -1.388] |
| Reward sensitivity | -0.108 | [-3.787, 4.526] |

The model-coded response mean refers to selection of option 2, not selection of
the higher-effort option.

## Ordinal and piecewise analyses

Overall ordinal accuracy is approximately 0.861 in MDD and 0.929 in controls.
The unadjusted difference is approximately -0.067; the age-, sex-, and success-
adjusted contrast is approximately -0.050.

The one-slope control model has BIC approximately 3103.10. In MDD, the best
two-segment candidate places the break after 30% MVC (BIC approximately 1700.34),
with only modest evidence relative to the one-slope alternative (Delta BIC 1.88).

