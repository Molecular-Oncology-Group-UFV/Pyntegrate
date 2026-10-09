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

# numpy 1.26 and pybedtools 0.9.1 only build on Python 3.9-3.12, so pick a
# compatible interpreter (override with PYTHON=/path/to/python).
if [ -z "$PYTHON" ]; then
    for candidate in python3.12 python3.11 python3.10 python3.9 python3; do
        if command -v "$candidate" >/dev/null 2>&1 &&
           "$candidate" -c 'import sys; sys.exit(not (3, 9) <= sys.version_info[:2] <= (3, 12))'; then
            PYTHON="$candidate"
            break
        fi
    done
fi
if [ -z "$PYTHON" ] ||
   ! "$PYTHON" -c 'import sys; sys.exit(not (3, 9) <= sys.version_info[:2] <= (3, 12))'; then
    echo "ERROR: Pyntegrate needs Python 3.9-3.12, but none was found." >&2
    echo "       Install one, or set PYTHON=/path/to/python3.10 and re-run." >&2
    exit 1
fi

echo "Creating virtual environment at $VENV_DIR with $("$PYTHON" --version) ..."
"$PYTHON" -m venv "$VENV_DIR"

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

pip install --upgrade pip
echo "Installing build tools (setuptools, wheel, numpy, Cython)..."
# pybedtools 0.9.1 cannot be built with setuptools >= 82; numpy/Cython match requirements.txt.
pip install "setuptools<82" wheel numpy==1.26.3 Cython==3.0.8

# pip's isolated build environments sometimes lack setuptools for legacy packages,
# so install pybedtools without build isolation to use the current venv's setuptools.
# (Not via conda: that would install it into the active conda env, not this venv.)
pip install --no-build-isolation pybedtools==0.9.1

# Finally install this package in editable mode
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
