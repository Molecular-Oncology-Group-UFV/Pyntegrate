# Installation

## Requirements

Pyntegrate works best inside an isolated Conda environment. Before installing, make sure you have:

| Requirement | Why |
|---|---|
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

=== ":material-rocket-launch: One-command script (recommended)"

    The repository ships a helper script that creates a `venv`, installs Pyntegrate in editable mode from a local clone, and installs `bedtools` automatically (via `conda`, `apt-get`, or `brew`, whichever is available) — handy for development or for running the [`examples/`](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples) scripts:

    ```bash
    git clone https://github.com/SalvaaGomez/Pyntegrate.git
    cd Pyntegrate
    ./setup_env.sh          # creates .venv/, installs Pyntegrate + bedtools
    # ./setup_env.sh /path/to/venv   # optional: custom venv location

    source .venv/bin/activate
    ```

=== ":material-package-variant: Conda (manual)"

    Creates a clean environment, installs the required runtime dependencies, and installs Pyntegrate from GitHub:

    ```bash
    conda create --name pyntegrateTest python=3.10 -y
    conda activate pyntegrate_final

    conda install -c conda-forge numpy=1.26.3 -y
    conda install -c bioconda bedtools -y
    conda install pandas -y
    conda install openpyxl -y
    conda install bioconda::samtools -y
    conda install conda-forge::pycirclize -y

    cd ~/pyntegrate
    git clone git@github.com:SalvaaGomez/Pyntegrate.git
    cd Pyntegrate/
    pip install -e .
    ```

=== ":material-cog: Plain venv (manual)"

    Uses Python's built-in `venv` instead of Conda. `bedtools` is not available via `pip`, so install it separately with your system package manager:

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate

    pip install --upgrade pip
    pip install git+https://github.com/SalvaaGomez/Pyntegrate.git

    sudo apt-get install bedtools  # or your platform's equivalent
    ```

## Install from a local copy

=== "From a ZIP download"

    ```bash
    conda create --name pyntegrateTest python=3.10 -y
    conda activate pyntegrateTest

    conda install -c conda-forge numpy=1.26.3 -y
    conda install -c bioconda bedtools -y

    unzip Pyntegrate-main.zip
    cd Pyntegrate-main
    pip install .
    ```

=== "From a Git clone"

    ```bash
    git clone https://github.com/SalvaaGomez/Pyntegrate.git
    cd Pyntegrate
    pip install .
    ```

!!! tip "Editable install for development"
    If you plan to modify Pyntegrate or run the [`examples/`](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples) scripts against a local checkout, install in editable mode instead: `pip install -e .`

## Optional dependencies

### HOMER

HOMER is optional and only needed for functions related to peak annotation and enrichment analysis.

Follow the [HOMER installation instructions](http://homer.ucsd.edu/homer/introduction/install.html) for your system.

## Verify the installation

After installation, check that the package imports correctly:

```python
import Pyntegrate
print(Pyntegrate.__version__)
```

If the import succeeds, the installation is ready. You can also confirm that `bedtools` is available:

```bash
which bedtools
bedtools --version
```

## Troubleshooting

??? failure "NumPy version conflicts"
    If you see an error related to NumPy, install or reinstall the version used by the project:

    ```bash
    conda install -c conda-forge numpy=1.26.3 -y
    ```

    Then reinstall Pyntegrate:

    ```bash
    pip install --force-reinstall git+https://github.com/SalvaaGomez/Pyntegrate.git
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
