# Installation

## Requirements

Pyntegrate works best inside an isolated Conda environment. Before installing, make sure you have:

| Requirement | Why |
|---|---|
| Python 3.9-3.12 | `numpy` 1.26 and `pybedtools` 0.9.1 do not build on Python 3.13+ (3.10 recommended) |
| `conda` | Manages the Python environment and native dependencies (e.g. `bedtools`) |
| Git | Needed to clone the repository or install directly from GitHub |
| C/C++ build toolchain | Required to compile some scientific-Python dependencies |
| `pip` | Installs Pyntegrate and its Python dependencies inside the environment |

=== "Linux"

    Install the basic system tools first:

    ```bash
    sudo apt update
    sudo apt install -y git build-essential
    ```

=== "macOS"

    Install the Apple command line tools:

    ```bash
    xcode-select --install
    ```

## Quick Installation

Pick the workflow that matches your setup. The setup script is the fastest path — it creates the environment, installs Pyntegrate, and takes care of `bedtools` for you in one command.

=== ":material-rocket-launch: One-command script"

    The repository ships a helper script that creates a `venv`, installs Pyntegrate in editable mode from a local clone, and installs `bedtools` automatically (via `conda`, `apt-get`, or `brew`, whichever is available) — handy for development or for running the [`examples/`](https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate/tree/main/examples) scripts:

    ```bash
    git clone https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate.git
    cd Pyntegrate
    ./setup_env.sh          # creates .venv/, installs Pyntegrate + bedtools
    # ./setup_env.sh /path/to/venv   # optional: custom venv location
    # PYTHON=/usr/bin/python3.11 ./setup_env.sh   # optional: choose the interpreter

    source .venv/bin/activate
    ```

=== ":material-package-variant: Conda (manual)"

    Creates a clean environment and installs Pyntegrate in editable mode from a local clone. Conda installs only what `pip` cannot (`bedtools`, `samtools`, and `pybedtools`, which often fails to compile under `pip`); `pip install -e .` then installs the remaining Python dependencies at the versions pinned in `requirements.txt`:

    ```bash
    conda create --name pyntegrate -c conda-forge -c bioconda \
        python=3.10 bedtools samtools pybedtools=0.9.1 -y
    conda activate pyntegrate

    git clone https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate.git
    cd Pyntegrate
    pip install -e .
    ```

    !!! warning "Activate the environment you just created"
        Make sure the name passed to `conda activate` matches the one used in `conda create`. Otherwise everything is installed into a different (possibly older) environment.


## Install from a local copy

=== "From a ZIP download"

    ```bash
    conda create --name pyntegrate -c conda-forge -c bioconda \
        python=3.10 bedtools samtools pybedtools=0.9.1 -y
    conda activate pyntegrate

    unzip Pyntegrate-main.zip
    cd Pyntegrate-main
    pip install .
    ```

!!! tip "Editable install for development"
    If you plan to modify Pyntegrate or run the [`examples/`](https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate/tree/main/examples) scripts against a local checkout, install in editable mode instead: `pip install -e .`

## Optional dependencies

### HOMER

HOMER is optional and only needed for functions related to peak annotation and enrichment analysis.

Follow the [HOMER installation instructions](http://homer.ucsd.edu/homer/introduction/install.html) for your system.

## Verify the installation

After installation, check from the terminal (with the environment activated) that the package imports correctly:

```bash
python -c "import Pyntegrate; print(Pyntegrate.__version__)"
```

If the import succeeds, the installation is ready. You can also confirm that `bedtools` is available:

```bash
which bedtools
bedtools --version
```

## Troubleshooting

??? failure "Conda solver errors (`excluded by strict repo priority`, unsatisfiable packages)"
    `conda-forge` must take priority over `bioconda`. Check the order with `conda config --show channels`; `conda-forge` should be listed first. To fix it:

    ```bash
    conda config --add channels bioconda
    conda config --add channels conda-forge   # the last one added goes on top
    ```

    Avoid `conda install`-ing packages that `pip install -e .` already installs (e.g. `numpy`, `pandas`, `pycirclize`): `pip` will replace them with the versions pinned in `requirements.txt`, and mixing both can cause conflicts.

??? failure "NumPy version conflicts"
    If you see an error related to NumPy, reinstall the version pinned by the project from inside the activated environment:

    ```bash
    pip install --force-reinstall numpy==1.26.3
    ```

??? failure "`bedtools` not found"
    If a function complains that `bedtools` is missing, verify that it is installed and available on your `PATH`:

    ```bash
    which bedtools
    ```

    If it is installed but not found, reactivate the Conda environment, or add the binary's location to `PATH` manually.

## Building the documentation locally

This site is built with [MkDocs](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/). Those tools are not part of Pyntegrate's runtime dependencies, so install them separately:

```bash
pip install -r docs/requirements.txt

mkdocs serve   # live preview at http://127.0.0.1:8000
mkdocs build   # static site in site/
```

---

Ready to go? Continue to [Getting Started](getting_started.md).
