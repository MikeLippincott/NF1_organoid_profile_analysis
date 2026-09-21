#!/bin/bash

git_root=$(git rev-parse --show-toplevel)

ENV_PATH="$git_root/.venv"


if [[ -e "$ENV_PATH" ]]; then
    echo "Using existing virtual environment at $ENV_PATH"
else
    mkdir -p "$ENV_PATH"
    echo "Created virtual environment directory at $ENV_PATH"
fi

# run twice to ensure we are in a clean environment
# and not accidentally using an existing one
# try but will not fail if not in a conda environment
loops=2
for _ in $(seq 1 $loops); do
    if command -v conda >/dev/null 2>&1; then
        conda deactivate >/dev/null 2>&1 || true
    fi
    deactivate >/dev/null 2>&1 || true
done
unset VIRTUAL_ENV

rm -f uv.lock
rm -rf .venv

uv venv
uv sync

# shellcheck disable=SC1091
source .venv/bin/activate
# Use RELATIVE path - simple and reliable
uv pip install -e ./utils

echo "Virtual environment setup complete. To activate it, run:"
echo "source $ENV_PATH/bin/activate"

# set up uvr
echo "Setting up uvr..."
# check if r version is installed and meets the requirement
# parse this from the uvr.toml file
R_VERSION_REQUIRED=$(grep -A 1 "r-version" uvr.toml | grep -o ">=.*")
ENV_NAME=$(grep -m1 -E '^name[[:space:]]*=' uvr.toml | sed -E 's/^name[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/')
R_VERSION_INSTALLED=$(R --version | head -n 1 | awk '{print $3}')
if [[ $(printf '%s\n' "$R_VERSION_REQUIRED" "$R_VERSION_INSTALLED" | sort -V | head -n1) != "$R_VERSION_REQUIRED" ]]; then
    echo "R version $R_VERSION_REQUIRED or higher is required. Installed version is $R_VERSION_INSTALLED."
    echo "Please install the required R version."
    exit 1
fi

if [[ ! -d ".uvr" ]]; then
    echo "Setting up uvr environment..."
else
    rm -rf .uvr
fi

uvr sync

echo "Setting up Jupyter kernel for the virtual environment..."

ENV_NAME=$(grep -m1 -E '^name[[:space:]]*=' uvr.toml | sed -E 's/^name[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/')

cat > /tmp/install_kernel.R << EOF
IRkernel::installspec(name = "$ENV_NAME", displayname = "$ENV_NAME")
EOF

uvr run /tmp/install_kernel.R
rm -f /tmp/install_kernel.R

KERNEL_DIR=$(jupyter kernelspec list | awk -v name="$ENV_NAME" 'tolower($1) == tolower(name) {print $2}')

if [[ -z "$KERNEL_DIR" ]]; then
    echo "Could not find installed kernelspec matching '$ENV_NAME'"
    jupyter kernelspec list
    exit 1
fi

PROJECT_LIB="$(pwd)/.uvr/library"
R_BIN="/home/lippincm/.uvr/r-versions/4.3.0/lib/R/bin/R"

cat > "$KERNEL_DIR/kernel-wrapper.sh" << EOF
#!/bin/bash
export R_LIBS="$PROJECT_LIB"
exec "$R_BIN" --slave -e 'IRkernel::main()' --args "\$1"
EOF
chmod +x "$KERNEL_DIR/kernel-wrapper.sh"

python3 - "$KERNEL_DIR/kernel.json" "$KERNEL_DIR/kernel-wrapper.sh" << 'PYEOF'
import json, sys
path, wrapper = sys.argv[1], sys.argv[2]
with open(path) as f:
    spec = json.load(f)
spec["argv"] = [wrapper, "{connection_file}"]
spec.pop("env", None)
with open(path, "w") as f:
    json.dump(spec, f, indent=2)
PYEOF

echo "Final kernel.json:"
cat "$KERNEL_DIR/kernel.json"
echo "Wrapper script:"
cat "$KERNEL_DIR/kernel-wrapper.sh"
