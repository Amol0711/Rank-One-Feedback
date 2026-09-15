Transient response simulations
==============================

This package evaluates a fixed three-state continuous-time feedback model in
raw and normalized quadratic state metrics. It contains numerical inputs,
response evaluations, stationary-time diagnostics, finite sampled gain
searches, directed interval bounds, and independent integer-arithmetic replay.
No external dataset is required.

The physical state order is theta, xi, i. The raw metric is the identity in
those physical coordinates. The normalized metric has diagonal entries
1, 1/25, 1/4. Both metrics remain fixed during each gain comparison.

Scope of reproduction
---------------------

The response workflow recalculates matrix-exponential samples at supplied
time grids. The stationary workflow evaluates all 22 gains in each metric
at both 80 and 110 decimal digits. The direction workflow evaluates five
unshifted and six shifted cases, transporting linked singular vectors in the
prescribed metric and finding the shifted stationary times independently.

The finite-gain workflow recalculates sampled searches, local candidates,
detail curves, spectra and high-precision stationary comparisons. These are
numerical candidates, not certified global gain optimizers. Time searches
and spectral crosschecks do not provide interval enclosures.

The interval workflow computes three fixed response comparisons at both
50 and 80 decimal interval digits. Each includes independently enclosed
time samples, time-zero coverage, an intersample factor, terminal contraction,
and every cell of a reciprocal-gain partition. The integer replay recomputes
the corresponding inequalities with a 256-bit fixed grid and degree-34
matrix-exponential expansion. Arithmetic implementations and their execution
remain part of the assurance assumptions.

Prediction coefficients, reference direction vectors, comparison multipliers,
and the comparison rates 34359738368 and 274877906944 are supplied inputs.
This package evaluates responses and residuals against those inputs. It does
not reconstruct the supplied comparison ranges or multipliers. The input
coefficients contain 220 decimal digits. Evaluations use at most 180 digits.

The supplied range observations and four-state reference-accuracy data are
retained numerical records and are exported without recomputation. The two
supplied gain sweeps and response-marker file are retained as numerical data;
the new stationary and finite-gain runs provide their own independent outputs.
The export command is a data export, not a solver execution. No single command
is claimed to regenerate every supplied observation from first principles.

Environment and quick start
---------------------------

The tested environment is Linux with Python 3.13.5. The exact dependency
versions are pinned in requirements.txt. Windows PowerShell and macOS/Linux
shell commands are provided below; the recorded numerical validation was
performed on Linux, not on every operating system.

Run these commands from the extracted package root. Create the environment
outside the package. After activation, python must refer to that environment.

macOS or Linux (Bash or Zsh):

    python3 -m venv ../response_environment
    source ../response_environment/bin/activate
    python -m pip install -r requirements.txt
    export PYTHONDONTWRITEBYTECODE=1
    export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    python -B verify_release.py
    python -B -m pytest -q -p no:cacheprovider tests
    python -B reproduce.py --mode check --output ../response_check
    python -B reproduce.py --mode responses --output ../response_samples

Windows PowerShell:

    py -3.13 -m venv ..\response_environment
    $PY = "..\response_environment\Scripts\python.exe"
    & $PY -m pip install -r requirements.txt
    $env:PYTHONDONTWRITEBYTECODE = "1"
    $env:PYTEST_DISABLE_PLUGIN_AUTOLOAD = "1"
    $env:OPENBLAS_NUM_THREADS = "1"
    $env:OMP_NUM_THREADS = "1"
    & $PY -B verify_release.py
    & $PY -B -m pytest -q -p no:cacheprovider tests
    & $PY -B reproduce.py --mode check --output ..\response_check
    & $PY -B reproduce.py --mode responses --output ..\response_samples

Dependency installation may use the network. Numerical execution itself does
not use the network or launch subprocesses. Python 3.13 must be installed for
the Windows command above. The direct interpreter path avoids any dependency
on a PowerShell activation-script execution policy.

Workflow selection
------------------

Use one mode at a time, or all for the complete sequence. From an activated
macOS/Linux environment:

    python -B reproduce.py --mode all --output ../response_results

In PowerShell after defining $PY and the environment variables above:

    & $PY -B reproduce.py --mode all --output ..\response_results

Mode             Operation
check            Verify stored records, scalar arithmetic, and schemas.
export           Copy retained CSV files and generate selected_diagnostics.csv.
responses        Recompute 1607 matrix-exponential response samples.
stationary       Evaluate 44 metric/gain cases at 80 and 110 decimal digits.
directions       Recompute five unshifted and six shifted direction cases.
comparison       Evaluate residuals at the supplied unshifted comparison rate.
finite-gain      Recompute sampled searches and high-precision local candidates.
intervals        Recompute three directed response comparisons at two precisions.
integer-replay   Recompute their inequalities with integer/fraction arithmetic.
all              Execute the eight modes from export through integer-replay.

The all mode also performs the common input checks. It can require substantially
more computation than check or responses. No estimated execution time is assumed.
Run python -B reproduce.py --help for the complete command syntax.

Use a NEW output directory for every run. Existing output directories, symlink
output paths, and outputs overlapping the input package are rejected. Rename
the requested output directory for a repeat run; never overwrite retained data.

Every run writes execution.json and file_access.json outside the package.
A successful run reports status "passed" and source_unchanged true. Each
calculation writes its numerical outputs and an execution log to its own output
subdirectory. Access records use PACKAGE, OUTPUT, DEPENDENCY and SYSTEM path
tokens rather than machine-specific absolute paths. A Python runtime audit hook
rejects network activity, subprocess creation and file access outside the
numerical package, output tree and Python installation. It is an execution
check, not an operating-system security boundary.

Inputs and outputs
------------------

configs/fixed_inputs.json defines the physical matrices, state metrics,
reference vectors, prediction coefficients and supplied comparison constants.
configs/certificate_cases.json specifies every rational incumbent, time step,
cutoff, witness seed, decay requirement and reciprocal partition.
configs/response_grids.json fixes the response sample times.
configs/diagnostic_selection.json is the single versioned selected-row policy.
It generates selected_diagnostics.csv without changing the full data grids.
The current view selects normalized cases. Its CSV output retains a metric
identifier for provenance; display_columns defines the six-column numeric
presentation. Both metrics and every original gain remain in the full data.
No second row-selection list is maintained in documentation.

data contains the complete numeric CSV products. reference contains retained
high-precision outputs and complete interval records. numerics contains only
fixed-model numerical implementations and response-bound arithmetic. tests
contains regression and adverse-input checks. DATA_DICTIONARY.md describes
all CSV fields and the structured numerical records. NUMERICAL_METHODS.txt
specifies numerical conventions and comparison tolerances.

The manifest covers all distributed files except the manifest itself.
verify_release.py rejects extra files, missing files and altered content.
A checkout's root .git entry is ignored as local version-control metadata;
it is not distributed or used by any numerical workflow. Python bytecode and
pytest caches are disabled by the commands above. Keep results and environments
outside the package.

Numerical comparison conventions
--------------------------------

Stationary heights, times, gain partials and signed residuals are compared
against retained high-precision values. Matching within 1e-55 is required for
the primary stationary and direction comparisons. Finite-gain comparisons
use 1e-50 for nonzero quantities and an absolute near-zero scale of 1e-60.
Floating-point response samples use a scaled tolerance of 2e-10. Directed
interval runs must recover the same strict displayed inequalities; identical
binary endpoints are not required across numerical backends.

Gain partials hold the evaluated time fixed and are obtained from block
matrix exponentials. A small height-relative discrepancy does not imply a
small correction-relative discrepancy. The correction quotient is omitted
when its denominator is at most 1e-8 times the reference level. Empty fields,
JSON null and the explicit plotting token nan denote omitted quantities,
not numerical zeros. Signed residuals are retained through cancellation.

Numerical-method identifiers
----------------------------

Matrix-exponential gain partials use DOI 10.1137/080716426. The numerical
interval-exponential context is recorded by the identifiers in
NUMERICAL_METHODS.txt. Dependencies are installed rather than vendored.

Repository integrity
--------------------

The distributed .gitattributes keeps text files in LF format after checkout.
The .gitignore prevents common local caches and environments from being staged,
but it does not exempt those files from numerical inventory verification. Keep
all environments and output directories outside this package.

README.md and README.txt contain the same instructions in rendered and plain-text
formats. Changing either one, adding a file, or changing a data value requires an
intentional update of MANIFEST.json. Never bypass a failed integrity check.
