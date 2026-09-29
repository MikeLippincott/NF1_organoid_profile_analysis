#!/bin/bash
# Launch the data_viewing Streamlit app from anywhere inside the repo.
# Usage: bash data_viewing/run_app.sh [port]   (default port: 8501)


PORT="${1:-8501}"
GIT_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

uv run streamlit run "$GIT_ROOT"/data_viewing/app.py --server.port "${PORT}" --theme.base light
