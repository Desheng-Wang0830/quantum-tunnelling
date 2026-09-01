# One-loop quantum tunnelling through finite delta arrays

This repository contains the MATLAB source code, archived numerical results,
plotting utilities, and validation material accompanying Desheng Wang's MSc
thesis on scalar-field tunnelling through finite arrays of planar delta layers.

The calculation retains the elastic correction through strict
\(O(\mu^2)\). The background scattering problem is treated exactly for a
finite delta array, while the loop contribution is assembled from the four
implemented placements \(G_1\)–\(G_4\) under the prescription stated in the
thesis.

## Start here

There are four distinct tasks. Choose the one matching what you want to check.

| Task | What is calculated? | Recommended starting point |
|---|---|---|
| Inspect the archived results | Nothing; open the supplied MAT, CSV, PNG, PDF, or FIG files | The four case archives |
| Regenerate tables and figures | Post-processing only; the loop integrals are not repeated | `replot_arbitrary_N_results.m` |
| Recompute an energy scan | The implemented four-placement numerical calculation is repeated | `run_arbitrary_N_energy_scan.m` with settings loaded from an archived MAT file |
| Run the frozen validation suite | The R00–R10 numerical and analytic tests are recomputed | `run_final_validation.m` |

The fastest independent check is to regenerate the CSV files and figures from
the archived MAT files. A complete recomputation of all four 200-point scans is
substantially more expensive.

## Repository contents

| Archive | Contents and purpose |
|---|---|
| `N_delta_1loop_correction.zip` | MATLAB driver 1.3.0, post-processing 1.3.0, numerical engine 1.9.6, validation suite, and the post-freeze heavy-mass driver |
| `N=1,g=0.3.zip` | Centred single-layer result, final MAT, restart checkpoint, CSV tables, and individual plots |
| `N=2,sys,g=0.3.zip` | Same-sign mirror-symmetric two-layer result: final MAT and restart checkpoint |
| `N=3,g=0.3.zip` | Same-sign mirror-symmetric three-layer result: final MAT and restart checkpoint |
| `N=2,asys,g=0.3,-0.3.zip` | Mixed-sign barrier–well dimer: final MAT and restart checkpoint |
| `chapter5data_set.zip` | Combined plotting scripts, rendered thesis figures, and the complete heavy-mass CSV/MAT/checkpoint record |

The historical filenames beginning with `chapter6_` were retained after the
numerical-results chapter was renumbered. They correspond to the results
figures in the current thesis.

Each case archive contains
`arbitrary_N_energy_G1_G4_pole_scan_results.mat`. This is the canonical final
result and contains a MATLAB variable named `result`, including the complete
settings structure used for that run. The much larger
`energy_scan_checkpoint.mat` is required only to resume an interrupted scan;
it is not needed to inspect or replot the final results.

## Archived production cases

All four primary scans use

- \(m_\phi=0.2\), \(m_\Phi=1.4\), and \(\mu=3\);
- 200 equally spaced external energies over \(0.24\leq E\leq1.10\);
- production mode with `orderFactor = 1.0`;
- quadrature orders
  `q4/qPerp/collective/cutRadial/cutAngle/endpoint = 20/18/24/20/20/20`;
- driver 1.3.0 and numerical engine 1.9.6.

| ID | Geometry | Positions | Strengths | Archived status |
|---|---:|---|---|---|
| C1 | \(N=1\) | `[0]` | `[0.3]` | 200/200 completed and accepted |
| C2 | \(N=2\), same sign | `[-6, 6]` | `[0.3, 0.3]` | 200/200 completed and accepted |
| C3 | \(N=3\), same sign | `[-6, 0, 6]` | `[0.3, 0.3, 0.3]` | 200/200 completed and accepted |
| C4 | \(N=2\), mixed sign | `[-6, 6]` | `[0.3, -0.3]` | 200/200 completed and accepted |

The value \(\mu=3\) is a display choice used to make the retained correction
visible in the thesis figures. Because the formulas are truncated at strict
\(O(\mu^2)\), algebraic \(\mu^2\) rescaling does not by itself bound omitted
\(O(\mu^4)\) terms or establish perturbative convergence.

For these explicit-geometry runs, the authoritative geometry is stored in
`result.profile.positions`, `result.profile.strengths`, and the corresponding
fields under `result.data`. The unused `geometry.totalStrength` field inherited
from the generic settings template should not be used to reconstruct an
explicit case.

## Requirements

The archived calculation was performed with

- MATLAB R2025b Update 5;
- Windows 64-bit (`PCWIN64` / `win64`);
- driver 1.3.0;
- engine 1.9.6;
- validation hardening 1.0.0.

The engine release contract declares MATLAB R2016b as its minimum version,
but the complete plotting and reproduction workflow documented here has been
checked only in MATLAB R2025b Update 5. The Parallel Computing Toolbox is
optional; the engine has a deterministic serial fallback. A standard desktop
MATLAB installation with a JVM is required for the SHA-256 source-manifest
checks.

Allow at least 1 GB of free space when extracting all four data archives.

## Recommended extraction layout

Download all six ZIP files, then arrange the extracted material as follows.
The names `chapter6_work` and `chapter6_figures` are only working-directory
names and may be changed.

```text
quantum-tunnelling/
├── N_delta_1loop_correction/
│   ├── run_arbitrary_N_energy_scan.m
│   ├── replot_arbitrary_N_results.m
│   ├── run_final_validation.m
│   ├── +nqscan/
│   ├── engine/
│   └── final_validation/
├── chapter6_work/
│   ├── N=1,g=0.3/
│   ├── N=2,sys,g=0.3/
│   ├── N=3,g=0.3/
│   └── N=2,asys,g=0.3,-0.3/
└── chapter6_figures/
    └── chapter6数据/
        ├── plot_chapter6_four_cases_tree_vs_strict_oneloop.m
        ├── plot_chapter6_four_cases_reflection_tree_vs_strict_oneloop.m
        ├── plot_chapter6_G2_G3_symmetry_comparison.m
        └── mPhi_decoupling_v130/
            ├── mPhi_decoupling_v130_data.csv
            ├── mPhi_decoupling_v130_results.mat
            └── mPhi_decoupling_v130_checkpoint.mat
```

Extract each of the four case ZIP files into `chapter6_work`. Each archive
already contains its own named case directory. Extract
`chapter5data_set.zip` into `chapter6_figures`.

## MATLAB setup

Start MATLAB and change directory to the inner code directory containing
`run_arbitrary_N_energy_scan.m`:

```matlab
codeRoot = 'C:\path\to\quantum-tunnelling\N_delta_1loop_correction';
cd(codeRoot)
clear functions
rehash
```

Do not begin with `addpath(genpath(...))`. The public entry points add the
required package and engine directories in a controlled order and verify that
MATLAB resolves the bundled engine rather than another copy on the path.

## 1. Regenerate one case without repeating the loop integrals

Set `caseRoot` to any one of the four extracted case directories:

```matlab
caseRoot = 'C:\path\to\quantum-tunnelling\chapter6_work\N=2,sys,g=0.3';

matPath = fullfile(caseRoot, ...
    'arbitrary_N_energy_G1_G4_pole_scan_results.mat');

[archivedResult, outputFiles] = replot_arbitrary_N_results( ...
    matPath, fullfile(caseRoot, 'results_replot_v1.3.0'));
```

This operation does **not** repeat a \(G_1\)–\(G_4\) integral. It reconstructs
the tables, source and term ledgers, singularity data, and figures from the
stored numerical arrays. The output is written to
`results_replot_v1.3.0`.

Important outputs include

```text
scan_summary.csv
barrier_profile.csv
G1_G4_term_contributions.csv
G1_G4_source_contributions.csv
G1_G4_Rose_four_source_contributions.csv
UV_unitarity_residuals_vs_energy.csv
q0_singularity_locations_vs_energy.csv
q0_branch_cut_snapshot_samples.csv
q0_branch_endpoints_and_collective_poles_vs_qperp.csv
collective_pole_contributions.csv
collective_pole_by_mode.csv
SCAN_REPORT.txt or REPLOT_REPORT_v1.3.0.txt
```

## 2. Rebuild the four-case tables and combined thesis figures

Only the C1 archive already contains the complete exported CSV collection.
The compact C2–C4 archives contain the final result MAT and checkpoint. First
generate the CSV files in all four case directories:

```matlab
repoRoot = 'C:\path\to\quantum-tunnelling';
codeRoot = fullfile(repoRoot, 'N_delta_1loop_correction');
dataRoot = fullfile(repoRoot, 'chapter6_work');

cd(codeRoot)
clear functions
rehash

caseNames = { ...
    'N=1,g=0.3', ...
    'N=2,sys,g=0.3', ...
    'N=3,g=0.3', ...
    'N=2,asys,g=0.3,-0.3'};

for k = 1:numel(caseNames)
    caseRoot = fullfile(dataRoot, caseNames{k});
    matPath = fullfile(caseRoot, ...
        'arbitrary_N_energy_G1_G4_pole_scan_results.mat');

    % Write the reconstructed CSV and figure files directly into the
    % extracted working copy of the case directory.
    replot_arbitrary_N_results(matPath, caseRoot);
end
```

Now add the directory containing the three combined plotting functions and
pass `dataRoot` explicitly:

```matlab
figureRoot = fullfile(repoRoot, ...
    'chapter6_figures', 'chapter6数据');
addpath(figureRoot, '-begin')

plot_chapter6_four_cases_tree_vs_strict_oneloop(dataRoot);
plot_chapter6_four_cases_reflection_tree_vs_strict_oneloop(dataRoot);
plot_chapter6_G2_G3_symmetry_comparison(dataRoot);
```

The scripts write PNG, vector PDF, and editable MATLAB FIG files next to the
plotting scripts. Passing `dataRoot` explicitly avoids relying on the older
internal expectation that the input archive is named `chapter6数据.zip`.

## Thesis-output map

| Current thesis output | Canonical archived input | Reproduction route |
|---|---|---|
| Four-case transmission comparison | C1–C4 `scan_summary.csv`: `T0Left`, `TLeftStrictOmu2` | `plot_chapter6_four_cases_tree_vs_strict_oneloop.m` |
| Four-case reflection comparison | C1–C4 `scan_summary.csv`: `R0Left`, `RLeftStrictOmu2` | `plot_chapter6_four_cases_reflection_tree_vs_strict_oneloop.m` |
| Symmetric \(G_1\)–\(G_4\) ledgers | C1–C3 `G1_G4_term_contributions.csv` | Replot C1–C3 from their final MAT files |
| Mixed-sign ledger and \(G_2/G_3\) comparison | C4 term ledger plus C2/C4 `scan_summary.csv` | Replot C4 and run `plot_chapter6_G2_G3_symmetry_comparison.m` |
| Heavy-mass decoupling figure | `mPhi_decoupling_v130_data.csv` and `mPhi_decoupling_v130_results.mat` | Inspect/replot the archived arrays, or rerun `run_mPhi_decoupling.m` with the recorded parameters |

### Heavy-mass archived result and full rerun

The heavy-mass directory contains the numerical source of Figure 5.5, not
only its rendered image. The MAT and CSV files record 26 completed masses
(25 logarithmically spaced values over \(0.5\leq m_\Phi\leq10^6\), plus the
inserted reference \(m_\Phi=1.4\)), signed \(\Delta T\) and \(\Delta R\), the
four term contributions, the fit mask, closure diagnostics, and fit
statistics. The run used

- explicit \(N=1\), position `[0]`, strength `[0.3]`;
- \(E=0.30\), \(m_\phi=0.2\), and \(\mu=3\);
- production `orderFactor = 1.0`;
- the 21 points with \(m_\Phi\geq5\) for the unweighted log-log fit.

The stored fit has slope `-1.9179634246402253`, regression standard error
`0.0119028948265120`, and \(R^2=0.999268759122075\). These are regression
statistics, not quadrature or tail-error bounds.

To inspect the archived result in MATLAB:

```matlab
heavyRoot = fullfile(figureRoot, 'mPhi_decoupling_v130');
payload = load(fullfile(heavyRoot, ...
    'mPhi_decoupling_v130_results.mat'), 'result');

result = payload.result;
loglog(result.mPhiValues, result.absoluteDeltaT, 'o-')
set(gca, 'XScale', 'log', 'YScale', 'log')
```

For a full recomputation, use the actual entry point
`run_mPhi_decoupling(codeRoot)`. Before running, set the editable base settings
to the explicit \(N=1\) configuration above, set `mPhiMaximum = 1e6` in the
heavy-mass driver, retain 25 logarithmic samples with the inserted reference
mass, and set `fitMinimumMPhi = 5`. The distributed defaults are deliberately
generic and therefore do not reproduce the thesis case without these edits.
The header comments inside this post-freeze driver retain the earlier working
name `run_mPhi_decoupling_left_v130`; the callable function and public filename
are `run_mPhi_decoupling` and `run_mPhi_decoupling.m`.

For provenance, use `result.profile.positions`, `result.profile.strengths`,
and the stored physics fingerprint. Some non-operative fields inherited by the
checkpoint request-key template retain the generic support or total-strength
summary; they were not the explicit profile passed to the engine for this run.

## 3. Recompute a complete archived energy scan

Calling `run_arbitrary_N_energy_scan()` with no argument does **not** reproduce
one of the four thesis cases. The default editable settings file is a generic
\(N=5\), \(\mu=0.1\) example.

For an exact case definition, load the settings stored in that case's final
MAT file. Change only the output directory and disable checkpoint reuse:

```matlab
repoRoot = 'C:\path\to\quantum-tunnelling';
codeRoot = fullfile(repoRoot, 'N_delta_1loop_correction');
caseMat = fullfile(repoRoot, 'chapter6_work', 'N=1,g=0.3', ...
    'arbitrary_N_energy_G1_G4_pole_scan_results.mat');

cd(codeRoot)
clear functions
rehash

archive = load(caseMat, 'result');
settings = archive.result.settings;

settings.execution.resume = false;
settings.output.directoryName = 'recomputed_C1_N1_g0p3';

freshResult = run_arbitrary_N_energy_scan(settings);
```

`settings.output.directoryName` must be a relative directory under the code
package. Use a different output directory for each case. The fresh computation
will write a new final MAT, checkpoint, tables, report, and requested figures.

The following comparison reports maximum absolute differences for several
primary stored arrays:

```matlab
referenceData = archive.result.data;
freshData = freshResult.data;

fields = { ...
    'externalEnergy', ...
    'treeProbabilities', ...
    'termDeltaProbabilities', ...
    'totalDeltaProbabilities'};

for k = 1:numel(fields)
    name = fields{k};
    difference = abs(referenceData.(name) - freshData.(name));
    fprintf('%-28s max abs difference = %.10e\n', ...
        name, max(difference(:)));
end
```

Cross-machine agreement should be assessed with declared absolute and relative
tolerances, not by requiring bitwise-identical MAT or PNG files.

### Direct inspection of the stored MATLAB data

The principal arrays can also be inspected without running the replotter:

| MATLAB field | Meaning |
|---|---|
| `result.data.externalEnergy` | 200-point external-energy grid |
| `result.data.treeProbabilities(:,1)` | Tree-level transmission \(T_0\) |
| `result.data.treeProbabilities(:,2)` | Tree-level reflection \(R_0\) |
| `result.data.totalDeltaProbabilities(:,1)` | Retained correction \(\Delta T\) |
| `result.data.totalDeltaProbabilities(:,2)` | Retained correction \(\Delta R\) |
| `result.data.correctedProbabilities` | Strict-order corrected transmission and reflection |
| `result.data.termDeltaProbabilities(:,j,c)` | \(G_j\) contribution for channel `c` |
| `result.data.sourceDeltaProbabilities` | Source-resolved probability ledger |
| `result.data.completed`, `accepted`, `numericalAccepted` | Per-energy completion and acceptance flags |

The four production archives are left-incidence production scans;
`result.data.fullSMatrixTested` is false for these stored runs. They should not
be described as amplitude-level full-S-matrix or perturbative-control scans.

## 4. Run quick code checks

These checks exercise the plotting/source grouping and auxiliary wrappers
without performing the complete production energy scans:

```matlab
cd(codeRoot)
clear functions
rehash

run_plot_wrapper_self_tests();
run_uv_audit_self_tests();
run_erman_feature_self_tests();
```

## 5. Run the frozen validation suite

Use the non-throwing entry point so that all reports are written even if a
numerical gate does not pass. Write the rerun output outside the frozen code
directory:

```matlab
repoRoot = 'C:\path\to\quantum-tunnelling';
codeRoot = fullfile(repoRoot, 'N_delta_1loop_correction');
validationOutput = fullfile(repoRoot, 'validation_rerun_output');

cd(codeRoot)
clear functions
rehash

validation = run_final_validation(validationOutput);
```

The saved reference run took approximately 6632 seconds (1 hour 50 minutes)
on its recorded machine. `run_final_validation_ci()` performs the same work but
throws an exception after writing the reports whenever the overall result is
not `PASS`; it is therefore not the recommended first run for this archive.

### Reference validation status

The reference report included in the code archive records

```text
Overall status: FAIL
9 PASS, 3 FAIL, 0 INCOMPLETE, 0 ERROR
```

The non-passing gates are:

| Test | Recorded result |
|---|---|
| R03, \(N=1\) analytic/numerical regression | The exact tree, \(\Gamma/\tau\), compact \(G_2\)–\(G_4\), probability, and collective-mode checks pass; the numerical \(G_1\) golden-kernel error is \(2.4668\times10^{-7}\), above the \(10^{-7}\) gate |
| R06, heavy-mass decoupling | The decreasing trend, fitted slopes, monotonicity, mapping-scale check, and loop closure pass; the largest masswise quadrature error is \(1.2677\times10^{-6}\), above the \(10^{-6}\) gate |
| R08, mapped-infinite versus finite-cutoff | The inner \(q_4\) comparison passes; the independent outer-\(q_\perp\) finite-cutoff comparison has not reached the required tail plateau |

The completion marker means that the suite finished and wrote its reports; it
does not mean that every gate passed.

### Source-manifest note for the public code ZIP

The frozen baseline manifest records 282 source and documentation files. All
270 baseline-listed files present in the public code ZIP match their recorded
hashes, but the ZIP omits the following 12 documentation files:

```text
ENGINE_PROVENANCE.txt
ERMAN_FEATURE_GUIDE_CN.md
FINAL_VALIDATION_CHECKLIST_CN.md
PACKAGE_VERIFICATION.txt
README_CN.md
RELEASE_NOTES_CN.md
RUN_FINAL_VALIDATION_CN.md
START_HERE_CN.txt
UV_UNITARITY_AUDIT_CN.md
final_validation/GAP_AUDIT_AND_COMPLETION_CN.md
final_validation/TOLERANCE_POLICY_CN.md
final_validation/reference_audit/INDEPENDENT_PYTHON_CROSSCHECK_CN.md
```

These are documentation files rather than executable numerical kernels. The
public ZIP also contains the post-freeze `run_mPhi_decoupling.m` driver, which
is not a baseline-listed file. Consequently, R99 reports a baseline-manifest
mismatch on this distribution even though the baseline-listed executable
files present still match their hashes. The decoupling driver and its archived
outputs are reproducibility additions rather than members of the frozen R99
source identity.

## Archive integrity

The SHA-256 digests of the files audited for this README are:

| File | SHA-256 |
|---|---|
| `N_delta_1loop_correction.zip` | `a3cc177e60e4fe3bb0fa957abf84321c5a2fd6f605e4e7d42d9ee8e7a05b9bab` |
| `N=1,g=0.3.zip` | `8888d63eb31b71499cfc86e85fbcbf1a362fab6231c05549d0aab42eee933ab8` |
| `N=2,sys,g=0.3.zip` | `15a257d92d3c0dc5e9b667e11346005fdcb425414c1542f4d1502528ea8190c5` |
| `N=3,g=0.3.zip` | `a3207b288a04cf13802cdea05a263d557e1acb4bf4e387387e2084352de6ec49` |
| `N=2,asys,g=0.3,-0.3.zip` | `8d52a35529d75d4fc69c9afcc2adc45a1decd4e9840310b4e3991a50c8977c2b` |
| `chapter5data_set.zip` | `8d36a23a75999da873d0114a329d79b7a5b7db396dfd4e5f9ad6a832b95628d4` |

The canonical final MAT files inside the four case archives have the following
digests:

| Case | Final result MAT SHA-256 |
|---|---|
| C1, \(N=1\) | `e967f0e5e3603e4580cd8232b9ef33160f99b8ee74771441fdc2307237a5a366` |
| C2, \(N=2\) same sign | `65df4d5dd8ba6cef997b77c2139b3aad77b762790e915a09050867eb68b51c3a` |
| C3, \(N=3\) same sign | `c4e522f5a12c64af0eb82e2c53f0ee89f9f3c572921a88c6ce76f207699f8f9a` |
| C4, \(N=2\) mixed sign | `c45c7b9e6105350a386a1ce8fa863f99422a5c218a27c85334520d3f5cabf30e` |

On Windows PowerShell, verify an archive with

```powershell
Get-FileHash '.\N_delta_1loop_correction.zip' -Algorithm SHA256
```

On Linux or macOS, use

```bash
sha256sum N_delta_1loop_correction.zip
```

## Known limitations and interpretation

1. **The heavy-mass thesis curve is reproducible, but its extreme tail is not
   a precision convergence claim.** `chapter5data_set.zip` supplies the exact
   mass grid, signed corrections, fit mask, term ledger, checkpoint, MAT/CSV
   arrays, and rendered figures. The supplied driver can recompute the scan
   after the recorded settings above are applied. Only one production order
   is archived, however, so the regression standard error does not bound
   quadrature, cancellation, or mapped-tail systematics at
   \(|\Delta T|\sim10^{-13}\). The independent R06 three-mass validation is a
   separate test rather than the source data for Figure 5.5.

   The archive also does not contain a single script that assembles every
   panel of the final composite ledger figures. Replotting produces the
   individual C1–C4 ledgers and the supplied comparison script generates the
   \(G_2/G_3\) panel; final multi-panel page assembly follows the thesis layout.

2. **The archived calculation is strict order \(O(\mu^2)\).** Reproduction of
   the stored arrays does not bound uncomputed \(O(\mu^4)\) contributions.

3. **The finite-array result is the implemented four-placement sum under the
   stated prescription.** The numerical archive should not be interpreted as
   an independent theorem establishing source completeness for every finite
   array, pole coalescence, or endpoint configuration.

4. **The stored historical environment files contain an absolute path such as
   `C:\Users\24701\...`.** This is provenance from the original Windows run,
   not a required runtime path. The executable public entry points determine
   their package paths at runtime.

5. **Raster-image hashes are not numerical acceptance tests.** Font rendering,
   graphics drivers, and MATLAB versions can change PNG bytes without changing
   the underlying arrays. Compare MAT/CSV numerical values instead.

## Troubleshooting

### MATLAB reports that the wrong engine is on the path

Return to the extracted code directory and reset cached functions:

```matlab
cd(codeRoot)
clear functions
rehash
```

Avoid adding a different N-delta package recursively to the same MATLAB path.

### A combined plotting script cannot find `scan_summary.csv`

The C2–C4 public archives contain compact MAT/checkpoint pairs. Run
`replot_arbitrary_N_results(matPath, caseRoot)` for each case first, writing
the reconstructed CSV files directly into the corresponding extracted working
case directory.

### A plotting script cannot find `chapter6数据.zip`

The public archive is named `chapter5data_set.zip`, while the older plotting
wrappers search for `chapter6数据.zip` in their no-input form. Extract the data
as described above and pass `dataRoot` explicitly to the plotting function.

### MATLAB reports that the output directory is invalid

The energy-scan driver requires `settings.output.directoryName` to be a
non-empty relative path under the code package. Absolute paths and paths
containing `.` or `..` are rejected.

## Citation

If these materials are used, please cite the accompanying MSc thesis and
record the exact Git commit or release tag together with the archive SHA-256
values above. The repository is maintained at
<https://github.com/Desheng-Wang0830/quantum-tunnelling>.
