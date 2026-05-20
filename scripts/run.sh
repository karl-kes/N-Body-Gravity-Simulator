#!/usr/bin/env bash
set -euo pipefail

# Simulation runner. Builds and runs main against tests/sim_ic.bin,
# which is produced by either jpl_compare.py fetch (solar) or
# generate_galaxies.py (galaxy).
#
# Usage:
#   ./scripts/run.sh                              # use existing sim_ic.bin
#   ./scripts/run.sh --fetch                      # fetch solar from JPL, then run
#   ./scripts/run.sh --galaxy --n 25000           # generate galaxy IC, then run
#   ./scripts/run.sh --force bh --theta 0.5       # pass through to main
#   ./scripts/run.sh --visualize                  # matplotlib viewer after sim
#   ./scripts/run.sh --render                     # Rerun dashboard after sim

FETCH=false
GALAXY=false
N=25000
START_DATE="1950-01-01"
MOONS="--moons"
FORCE="direct"
THETA="0.5"
VISUALIZE=false
RENDER=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --fetch)        FETCH=true;       shift ;;
        --galaxy)       GALAXY=true;      shift ;;
        --n)            N="$2";           shift 2 ;;
        --start)        START_DATE="$2";  shift 2 ;;
        --no-moons)     MOONS="";         shift ;;
        --force)        FORCE="$2";       shift 2 ;;
        --theta)        THETA="$2";       shift 2 ;;
        --visualize)    VISUALIZE=true;   shift ;;
        --render)       RENDER=true;      shift ;;
        -h|--help)
            echo "Usage: ./scripts/run.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --fetch            Fetch solar IC from JPL Horizons before running"
            echo "  --galaxy           Generate galaxy IC before running"
            echo "  --n VALUE          Particles per galaxy (default: 25000)"
            echo "  --start DATE       Solar fetch start date (default: 1950-01-01)"
            echo "  --no-moons         Solar: planets only (10 bodies)"
            echo "  --force MODE       direct or bh (default: direct)"
            echo "  --theta T          BH opening angle (default: 0.5)"
            echo "  --visualize        Open Matplotlib viewer after simulation"
            echo "  --render           Open Rerun dashboard after simulation"
            echo "  -h, --help         Show this help"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Stage IC if requested.
if [ "$FETCH" = true ] && [ "$GALAXY" = true ]; then
    echo "ERROR: cannot combine --fetch and --galaxy"
    exit 1
fi

if [ "$FETCH" = true ]; then
    echo "Fetching solar IC..."
    python3 src/jpl_compare.py fetch --start "$START_DATE" $MOONS
elif [ "$GALAXY" = true ]; then
    echo "Generating galaxy IC..."
    python3 src/generate_galaxies.py --n "$N"
fi

# Build:
echo ""
echo "Building..."
cmake -B build -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
cmake --build build --config Release

# Run:
echo ""
echo "Running simulation (force=$FORCE, theta=$THETA)..."
./build/main --force "$FORCE" --theta "$THETA"

# Validate (solar only):
if [ "$GALAXY" = false ]; then
    echo ""
    echo "Validating against JPL Horizons..."
    python3 src/jpl_compare.py compare || true
fi

# Visualize:
if [ "$VISUALIZE" = true ]; then
    echo ""
    if [ "$GALAXY" = true ]; then
        python3 src/visualize_galaxies.py
    else
        python3 src/visualize.py
    fi
fi

if [ "$RENDER" = true ]; then
    echo ""
    python3 src/render.py
fi

echo ""
echo "Done."