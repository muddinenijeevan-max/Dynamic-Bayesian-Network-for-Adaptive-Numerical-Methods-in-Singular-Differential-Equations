# Contributing

1. Fork the repository and create a feature branch.
2. Keep numerical solver code in `hofid_bvp.py`, benchmark definitions and output handling in `run_all_problems.py`, and plotting in `plot_all_problems.py`.
3. Install dependencies with `python -m pip install -r requirements.txt`.
4. Run the built-in solver checks with `python hofid_bvp.py` and test a small benchmark before opening a pull request.
5. Describe changes to numerical methods, convergence behavior, output formats, or dependencies in the pull request.
6. Do not commit generated result files, virtual environments, datasets, or machine-specific paths.
