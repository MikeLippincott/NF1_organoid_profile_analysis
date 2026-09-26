#!/bin/bash
# Launch the data_viewing Streamlit app from anywhere inside the repo.
# Usage: bash data_viewing/run_app.sh [port]   (default port: 8501)

# set -euo pipefail

PORT="${1:-8501}"
GIT_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

if [ -x "${GIT_ROOT}/.venv/bin/streamlit" ]; then
    STREAMLIT="${GIT_ROOT}/.venv/bin/streamlit"
else
    echo "streamlit not found in ${GIT_ROOT}/.venv. Run: source uv_setup.sh" >&2
    exit 1
fi

cd "${GIT_ROOT}" || exit 1
exec "${STREAMLIT}" run data_viewing/app.py --server.port "${PORT}" --theme.base light
