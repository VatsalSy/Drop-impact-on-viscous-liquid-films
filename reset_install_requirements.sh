#!/bin/zsh
# Build the Basilisk source committed alongside this project.
set -eu

if [[ "${1:-}" == "--hard" ]]; then
    hard=true
elif [[ $# == 0 ]]; then
    hard=false
else
    printf 'Usage: %s [--hard]\n' "$0" >&2
    exit 2
fi

root=$(cd -- "$(dirname -- "$0")" && pwd -P)
basilisk="$root/basilisk/src"
if [[ ! -f "$basilisk/Makefile" || ! -f "$basilisk/qcc.c" ]]; then
    printf '%s\n' 'Restore the vendored basilisk/ directory from this project checkout before building.' >&2
    exit 1
fi

(
    cd -- "$basilisk"
    if [[ "$OSTYPE" == darwin* ]]; then
        platform_config=config.osx
    else
        platform_config=config.gcc
    fi
    # The checkout carries the config link and binaries of the machine that
    # committed them. Keep a custom config; point a stock one at this platform.
    current_config=$(readlink config 2>/dev/null || true)
    if [[ ! -e config && ! -L config ]] ||
       [[ "$current_config" == config.osx || "$current_config" == config.gcc ]]; then
        ln -sfn "$platform_config" config
    fi
    # Rebuild from source for --hard, or when the committed binaries cannot run
    # here. "make -B" rebuilds ast/interpreter from the wrong directory and fails.
    if [[ "$hard" == true ]] || ! ./qcc --version >/dev/null 2>&1; then
        make clean
        rm -f ast/grammar  # not removed by "make clean"; make rebuilds it
    fi
    make
)

"$basilisk/qcc" --version
# Paths must remain literal even when the project directory contains spaces.
{
    printf 'export BASILISK=%q\n' "$basilisk"
    printf '%s\n' 'export PATH="$BASILISK:$PATH"'
} > "$root/.project_config"
