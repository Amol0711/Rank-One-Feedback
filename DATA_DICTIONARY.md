# Numerical data dictionary

The data use fixed physical coordinates and prescribed metrics. Null, an empty CSV field, and the plotting token `nan` indicate an omitted quantity. They never denote zero. The full data grids remain intact when a smaller diagnostic selection is exported.

## CSV inventory

| File | Rows | Numerical role |
|---|---:|---|
| `asymptotic_comparison.csv` | 44 | All 44 metric/gain rows at 110 digits. Recomputed by stationary. |
| `control_detail.csv` | 179 | 179 numerical detail-curve samples. Recomputed within finite-gain. |
| `control_spectrum.csv` | 385 | 385 direct spectral samples. Recomputed within finite-gain. |
| `direction_diagnostics.csv` | 5 | Five unshifted linked-direction cases. Recomputed by directions. |
| `peak_diagnostics.csv` | 22 | Normalized discrepancy subset of the same grid. |
| `response_k10.csv` | 402 | Fixed response-time grid, recomputed by responses. |
| `response_k100.csv` | 401 | Fixed response-time grid, recomputed by responses. |
| `response_k1000.csv` | 402 | Fixed response-time grid, recomputed by responses. |
| `response_k50.csv` | 402 | Fixed response-time grid, recomputed by responses. |
| `response_peak_markers.csv` | 4 | Four supplied stationary marker values. Exported; stationary produces fresh matching gain cases. |
| `supplied_range_observations.csv` | 6 | Six supplied scalar observations. Exported, not recomputed. |
| `supplied_reference_accuracy.csv` | 3 | Three supplied four-state accuracy rows. Exported, not recomputed. |
| `sync_gain_sweep.csv` | 96 | Supplied gain sweep. Exported without a new search. |
| `sync_peak_layer.csv` | 34 | Supplied large-gain sweep. Exported without a new search. |
| `weighted_envelope_diagnostics.csv` | 6 | Six shifted cases at their own times. Recomputed by directions. |

## Field definitions

| Field | Meaning |
|---|---|
| `EG_percent` | Supplied height discrepancy percentage using the reference height denominator. |
| `Et_percent` | Supplied time discrepancy percentage using the reference-time denominator. |
| `G` | Supplied reference height in the four-state data, not recomputed by the public workflows. |
| `Ghat` | Supplied predicted height in the four-state data. |
| `L` | Logarithmic prediction input log(chi*r/a_p). |
| `P` | Metric label of the supplied numerical reference, I4. |
| `R_D` | Signed sensitivity residual multiplied by r^3/(d*L^2); omitted if L = 0. |
| `R_G` | Signed height residual multiplied by r^2/L^2; omitted if L = 0. |
| `R_t` | Signed time residual multiplied by r^2/L; omitted if L = 0. |
| `a_iso` | Fixed unit decay parameter of the supplied numerical reference. |
| `a_p` | Supplied positive numerical loss coefficient for the fixed metric. |
| `alpha` | Largest real part of the eigenvalues of the fixed closed-loop matrix. |
| `assurance` | Numerical evidence classification and limitations attached to the row. |
| `below_supplied_comparison_range` | For normalized cases, whether the gain is below the supplied rate. Null denotes an unassigned range flag for the raw metric. |
| `beta` | Numeric margin input attached to the supplied scalar observation. |
| `bracket_lower` | Lower stationary-search bracket in fast time. |
| `bracket_upper` | Upper stationary-search bracket in fast time. |
| `case` | Identifier of a supplied scalar numerical case. |
| `chi` | Supplied positive fast-response coefficient. |
| `comparison_112L_over_r` | Comparison scale evaluated using the supplied multiplier 112. |
| `comparison_to_error_ratio` | Supplied direction comparison scale divided by the observed direction error. |
| `correction_defined` | True only when the correction denominator exceeds 1e-8*gamma. |
| `correction_denominator` | Absolute value of gamma minus predicted_peak. |
| `correction_error_pct` | 100*abs(peak-predicted_peak)/abs(gamma-predicted_peak), subject to the stated omission rule. |
| `correction_ratio_to_gamma` | Correction denominator divided by gamma. |
| `curvature` | Second time derivative at the stationary point. Negative values indicate a strict local maximum. |
| `d` | Positive rate-to-gain scale of the fixed numerical family. |
| `delta` | Nonnegative exponential weighting rate; zero or one in the retained direction data. |
| `denominator_policy` | Explicit reference denominators used for the supplied accuracy percentages. |
| `derivative_relative_defined` | True when the numerical gain partial is nonzero. |
| `dimension` | State dimension of the supplied numerical reference. |
| `direction_error` | Sum of the Euclidean errors of the linked pair in metric coordinates, with one simultaneous sign choice. |
| `dps` | Decimal working precision of the numerical calculation. |
| `evaluation_digits` | Working precision recorded for the supplied observation. |
| `fast_time` | Gain multiplied by model time for this fixed d = 1 family. |
| `gain` | Nonnegative scalar feedback gain used for this numerical row. |
| `gain_derivative` | Fixed-time gain partial evaluated at the stationary numerical time. |
| `gamma` | Supplied reference amplification level for the fixed metric. |
| `gap_below_gamma` | gamma minus stationary_height. |
| `h_p` | Supplied drift-dependent height-correction coefficient. |
| `height_error_pct` | 100*abs(peak-predicted_peak)/peak. |
| `height_residual` | Signed peak minus predicted_peak. |
| `in_supplied_comparison_range` | Whether this gain is at least its supplied comparison rate. |
| `j_p` | Supplied fixed numerical drift coefficient. |
| `log_parameter_below_one` | Boolean flag for L < 1. It is a numerical indicator, not a newly computed admissibility range. |
| `logarithmic_time` | Logarithmic prediction input for the selected weighting rate. |
| `lower_log2` | Supplied lower endpoint of a base-two logarithmic rate bracket. |
| `lower_rate` | Supplied lower rate endpoint; this bracket is not recomputed or outward-certified here. |
| `metric` | Fixed performance metric. raw denotes identity in physical coordinates and normalized denotes diagonal (1,1/25,1/4). |
| `norm_policy` | Label of the norm convention of the supplied scalar observation. No constructor is included. |
| `norm_witness` | Response norm evaluated at the specified time, without global-maximality certification. |
| `peak` | Stationary response singular value, or sampled/refined response candidate as specified by the file. |
| `peak_time` | Positive stationary time corresponding to peak. |
| `peak_witness` | Sampled and refined response-height candidate. |
| `predicted_gain_derivative` | Sensitivity prediction evaluated from supplied coefficients. |
| `predicted_peak` | Height prediction evaluated from the supplied coefficients and gain. |
| `predicted_time` | Time prediction L/r evaluated from supplied coefficients. |
| `projection_error_pct` | 100*abs(peak-gamma)/peak. |
| `r` | Rate d times gain. The fixed model uses d = 1. |
| `rate` | Rate at which the supplied four-state reference was evaluated. |
| `scaled_direction_error` | Observed direction error multiplied by r/L. |
| `sensitivity_error_pct` | 100*abs(gain_derivative-predicted_gain_derivative)/abs(gain_derivative), omitted at a zero denominator. |
| `sensitivity_opposite_sign` | Whether the numerical and predicted fixed-time gain sensitivities have opposite signs. |
| `sensitivity_residual` | Signed gain_derivative minus predicted_gain_derivative. |
| `singular_gap` | Leading singular value minus the next singular value. |
| `spectral_propagator_relative_error` | Relative difference between the direct and spectral matrix exponential. |
| `spectral_sensitivity_relative_error` | Relative difference between block-exponential and spectral-divided-difference gain partials. |
| `stationarity_residual` | Maximum of stationarity and linked singular-pair residuals in the numerical calculation. |
| `stationary_height` | Response height at its own stationary time, including the stated exponential weighting. |
| `stationary_time` | Stationary time of the unshifted or weighted response, as indicated by delta. |
| `status` | Evidence classification of a row. Stationary numerical witnesses are not interval-certified global extrema. |
| `supplied_comparison_rate` | Declared comparison-rate input. It is not recomputed by this package. |
| `t` | Supplied reference time in the four-state data. |
| `that` | Supplied predicted time in the four-state data. |
| `time` | Model-time coordinate, with its numerical role specified by the file. |
| `time_error_pct` | 100*abs(peak_time-predicted_time)/abs(peak_time). |
| `time_residual` | Signed peak_time minus predicted_time. |
| `time_witness` | Supplied sampled/refined candidate time. It is not an interval enclosure of a global maximizing time. |
| `upper_log2` | Supplied upper endpoint of a base-two logarithmic rate bracket. |
| `upper_rate` | Supplied upper rate endpoint; this bracket is not recomputed or outward-certified here. |

## Complete structured records

`reference/stationary.json` retains both precision runs, their row-level numerical quantities, local residual checks and matched-precision differences. Its `constants` member is numerical input/output data, not a general coefficient-construction routine.

`reference/directions.json` retains both precision runs for all eleven cases, the linked `u` and `v` vectors, their transported physical initial/response vectors, norm and sign checks, and the independent unshifted-time controls. `checks` entries ending in `residual` quantify numerical identities. Curvature, singular-gap and before/after derivative fields retain their local-stationarity interpretation.

`reference/comparison_residual.json` retains the two numerical evaluations at the supplied unshifted comparison rate. Each `bound_comparisons` member records the absolute observed residual, the comparison scale evaluated from the supplied multiplier, and their ratio. `correction_error_pct` is null in this case. No residual interval enclosure is attached.

`reference/finite_gain.json` retains the two sampled search densities, endpoint and interior candidates, 80/110-digit stationary rows, direct spectral checks and candidate curvature diagnostics. Conditional local curvature and stationary residuals do not certify a global minimizing gain.

`reference/interval_50/` and `reference/interval_80/` retain each complete directed numerical result. The `upper_certificate` contains gain, metric, time step, every sample, the grid maximum including time zero, the intersample numerical-abscissa bound and terminal contraction. The `exclusion_certificate` contains cutoff, witness-time parameter, rational seed vectors, normalization convention and all closed reciprocal-gain cells. Every interval has decimal `lower`/`upper` strings and exact binary endpoint tuples. `strict_margin` records rational evaluation of the fixed cubic decay conditions.

Binary endpoint tuples have sign, integer mantissa, binary exponent and bit count. Decimal endpoints are serialized outwards from the exact binary values. Small differences between fresh proposal endpoints do not invalidate a comparison when the directed strict inequalities are recovered.

## Generated selections and execution records

`configs/diagnostic_selection.json` alone defines the selected-row identities, order, source and columns. Export joins those identities to the full grid and the separate comparison row. It writes `selected_diagnostics.csv` and preserves all percentage definitions, signed residuals and omission flags. The version-2 policy retains metric identifiers in the selected CSV and specifies six numeric presentation columns through `display_columns`. Both complete metric grids remain unchanged. Documentation does not maintain a competing selection list.

`execution.json` distinguishes stored-record checking, data export, stationary numerical evaluation, sampled gain search and directed numerical verification. `file_access.json` uses path tokens rather than host-specific absolute names. Timing fields describe the local execution and are not reference tolerances.
