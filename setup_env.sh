#!/bin/bash
#
# Create a virtual environment, install Pyntegrate into it (editable mode),
# and make sure bedtools is available.
#
# Usage:
#   ./setup_env.sh [venv_dir]
#
# venv_dir defaults to .venv in the repo root.

set -e

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${1:-$HERE/.venv}"

echo "Creating virtual environment at $VENV_DIR ..."
python3 -m venv "$VENV_DIR"

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

pip install --upgrade pip
pip install "setuptools<82" wheel numpy Cython
pip install --no-build-isolation pybedtools==0.9.1
pip install -e "$HERE"

if command -v bedtools >/dev/null 2>&1; then
    echo "bedtools already available: $(command -v bedtools)"
elif command -v conda >/dev/null 2>&1; then
    echo "Installing bedtools via conda ..."
    conda install -c bioconda bedtools -y
elif command -v apt-get >/dev/null 2>&1; then
    echo "Installing bedtools via apt-get (may prompt for sudo password) ..."
    sudo apt-get update && sudo apt-get install -y bedtools
elif command -v brew >/dev/null 2>&1; then
    echo "Installing bedtools via Homebrew ..."
    brew install bedtools
else
    echo "WARNING: could not find conda, apt-get, or brew to install bedtools automatically."
    echo "         Install it manually, e.g.: conda install -c bioconda bedtools"
fi

cat <<EOF

Done. Pyntegrate is installed in editable mode in:
  $VENV_DIR

Activate the environment with:
  source $VENV_DIR/bin/activate
EOF
