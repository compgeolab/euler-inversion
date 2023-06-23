# Reproducing our results

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

```
git clone https://github.com/compgeolab/NAME_OF_THIS_REPOSITORY
```

or downloading a zip archive from the provided DOI in the `README.md`.

All source code used to generate the results and figures in the paper are in
the `code` folder. There you can find the [Jupyter](https://jupyter.org/)
notebooks and `.py` modules that perform the calculations and generate the
figures and results presented in the paper.

The LaTeX sources for the manuscript text and figures are in `paper`.

## Setting up your environment

You'll need a working Python 3 environment with all the standard
scientific packages installed (numpy, pandas, scipy, matplotlib, etc).
The easiest (and recommended) way to get this is to download and install the
[Mambaforge Python distribution](https://github.com/conda-forge/miniforge#mambaforge).
This will give you access to a Python 3 environment and the excellent `mamba`
package manager (a drop-in replacement for `conda`).

Besides the standard scientific packages, you'll also need to install some
extra libraries like: Numba for just-in-time compilation; Harmonica, Verde,
Boule and Pooch from the [Fatiando a Terra](https://www.fatiando.org) project;
PyGMT for generating maps and more.

Instead of manually installing them, they can all be automatically installed
using a conda environment.

1. Inside the cloned repository (or an unzipped version), create a new virtual
   environment from the `environment.yml` file by running:
   ```
   mamba env create -f environment.yml
   ```
1. Check the environment name in the `environment.yml` file.
1. Activate the new environment by running:
   ```
   mamba activate ENVIRONMENT_NAME
   ```

## Generating the results from the paper

All results and their associated figures are created by the notebooks in the
`code` folder. To fully reproduce them, run each notebook from top to bottom
following the numerical order of notebook file names.

With the environment activated, start a JupyterLab server by running:

```
jupyter lab
```

Once in JupyterLab, open each notebook in turn and run them from top to bottom.

> Some notebooks might take several minutes to run, depending on the resources
> of your system.

## Generating the preprint PDF

The `paper` folder provides a `Makefile` that includes commands that
automatically build the preprint PDF from the LaTeX source files using
[Tectonic](https://tectonic-typesetting.github.io/). GNU Make should be
installed in your conda environment from above.

Inside the `paper` folder, run the `make` command to generate `preprint.pdf`.
