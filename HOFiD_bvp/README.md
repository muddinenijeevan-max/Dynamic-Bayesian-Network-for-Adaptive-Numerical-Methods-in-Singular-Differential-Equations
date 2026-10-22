# HOFiD_bvp — Adaptive High-Order Finite-Difference Solver

A Python implementation of a high-order finite-difference method for singular and nonlinear two-point boundary-value problems (BVPs), with adaptive mesh refinement, damped Newton iteration, numerical diagnostics, saved results, and plotting utilities.

## Features

- High-order finite-difference approximations for first and second derivatives.
- Sparse operator matrices and sparse LU solves.
- Damped Newton iteration for nonlinear BVPs.
- Adaptive mesh refinement with block grading.
- Error estimation using solutions of order `p` and `p + 2` on the same mesh.
- Four benchmark problem families, including a problem with a known exact solution.
- Interactive selection of a benchmark problem and initial mesh-point count.
- CSV, JSON, and NPZ output, plus PNG plots.
- Reuses saved solutions when the matching result files are available.

## Requirements

- Python 3.10 or newer (Python 3.10+ recommended).
- NumPy
- SciPy
- Matplotlib
- pandas

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Project layout

```text
HOFiD_bvp/
├── hofid_bvp.py           # Numerical solver and diagnostics
├── run_all_problems.py    # Problem definitions, interactive runner, output
├── plot_all_problems.py   # Plot generation; reuses saved results
├── requirements.txt
├── .gitignore
├── LICENSE
├── CITATION.cff
├── README.md
└── results/               # Generated locally; ignored by Git
```

## Quick start

Run commands from the project root.

### 1. Solve benchmark problems

```bash
python run_all_problems.py
```

The script asks which problem to solve (`1`, `2`, `3`, `4`, or all, depending on the prompt) and whether to override the default initial mesh-point count. It then runs the selected experiments and saves their outputs under `results/`.

### 2. Plot saved results

```bash
python plot_all_problems.py
```

The plotting script asks for the problem selection and initial mesh setting. It loads matching saved results where possible. If a matching result is missing, it invokes the runner to compute and save it before plotting.

### 3. Run solver validation tests

```bash
python hofid_bvp.py
```

This runs the solver's built-in validation/self-tests, if enabled by the script.

## Generated outputs

Each run uses a folder named from the selection, for example:

```text
results/
└── problem_1_npts_15/
    ├── inputs.json
    ├── tables/
    │   ├── <experiment>_initial.csv
    │   └── <experiment>_final.csv
    ├── data/
    │   ├── <experiment>.npz
    │   └── <experiment>_summary.json
    ├── diagnostics/
    │   └── <experiment>_coefficients.csv
    └── figures/
        └── <generated plots>.png
```

Exact files depend on the selected problem and experiments. Generated output is intentionally excluded from version control so the repository stays lightweight. The empty `results/` directory is retained using `.gitkeep`.

## Numerical notes

- Mesh-point counts include both endpoints; intervals equal points minus one.
- The runner validates the requested initial mesh size against the selected experiment's minimum and maximum.
- The adaptive solver estimates error by comparing the order-`p` and order-`p + 2` solutions on the same mesh.
- Actual error against an exact solution is available only for benchmark cases that define one.
- Numerical convergence and accuracy depend on the problem, initial guess, order, tolerance, and mesh.

## Reproducibility

For reproducible experiments, record:
- Python and dependency versions;
- selected problem(s);
- initial mesh-point count (or default);
- solver order and tolerance;
- convergence status and output summary JSON.

Results are generated locally and may differ slightly across environments or dependency versions.

## Troubleshooting

**`ModuleNotFoundError` for NumPy, SciPy, Matplotlib, or pandas**

```bash
python -m pip install -r requirements.txt
```

**The wrong Python environment is being used**

```bash
python --version
python -m pip --version
```

Use the same Python executable for installing packages and running the scripts. In VS Code, select the interpreter for the environment where dependencies are installed.

**Plots are not generated**

Check that the selected experiment converged and inspect the terminal output and `results/<run>/data/` summaries.

## Citation

If you use this code in academic work, cite the repository and any source paper or method that the implementation is based on. Update `CITATION.cff` with the final author, repository URL, and publication details before release.

## License

This project is distributed under the MIT License; see [`LICENSE`](LICENSE).
