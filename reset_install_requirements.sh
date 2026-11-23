#!/bin/zsh
# Build the Basilisk source committed alongside this project.
set -eu

if [[ "${1:-}" == "--hard" ]]; then
    make_flags=(-B)
elif [[ $# == 0 ]]; then
    make_flags=()
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
    if [[ ! -e config && ! -L config ]]; then
        if [[ "$OSTYPE" == darwin* ]]; then
            ln -s config.osx config
        else
            ln -s config.gcc config
        fi
    fi
    make "${make_flags[@]}"
)

"$basilisk/qcc" --version
# Paths must remain literal even when the project directory contains spaces.
{
    printf 'export BASILISK=%q\n' "$basilisk"
    printf '%s\n' 'export PATH="$BASILISK:$PATH"'
} > "$root/.project_config"
