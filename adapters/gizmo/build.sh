#!/bin/bash
# cofactor adapter for GIZMO: an isolated build copy with a test-only hook, Config.sh from the given flags, no-FMA make.
# usage: build.sh SRC BUILD HOOK|- MAKE_JOBS "OPTIMIZE" -- CONFIG_FLAG...
set -eu
src=$1 build=$2 hook=$3 jobs=$4 opt=$5
shift 5
[ "${1:-}" = "--" ] && shift
mkdir -p "$build"
rsync -a --exclude .git --exclude 'output*' --exclude '*.o' --exclude GIZMO --exclude 'test/*/run_*' "$src"/ "$build"/
cd "$build"
if [ "$hook" != "-" ]; then
    if patch -p1 --dry-run -R -s -f < "$hook" > /dev/null 2>&1; then
        echo "cofactor: hook already present in the source tree"
    else
        patch -p1 -s --forward < "$hook" || { echo "cofactor: hook patch does not apply to this source tree"; exit 5; }
    fi
fi
printf '%s\n' "$@" > Config.sh
echo "== Config.sh"; cat Config.sh
make clean > /dev/null 2>&1 || true
make -j "$jobs" OPTIMIZE="$opt"
test -x GIZMO
