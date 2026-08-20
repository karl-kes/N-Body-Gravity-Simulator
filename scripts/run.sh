#!/usr/bin/env bash
set -euo pipefail

# Simulation Pipeline Runner
# Generates the IC, builds, runs, validates, and visualizes.
#
# Usage:
#   ./run.sh                                # use existing tests/sim_ic.bin
#   ./run.sh --fetch                        # fetch solar IC from JPL, then run
#   ./run.sh --galaxy                       # generate MW-Andromeda IC, then run
#   ./run.sh --galaxy --n 50000             # 50k particles per galaxy
#   ./run.sh --force bh --theta 0.5         # Barnes-Hut force at theta = 0.5
#   ./run.sh --galaxy --force bh --visualize
#   ./run.sh --fetch --visualize --render   # full solar pipeline with both viewers

FETCH=false
GALAXY=false
N=25000
START="1950-01-01"
NO_MOONS=false
FORCE="direct"
THETA=0.5
VISUALIZE=false
RENDER=false
FP32=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --fetch)        FETCH=true;       shift ;;
        --galaxy)       GALAXY=true;      shift ;;
        --n)            N="$2";           shift 2 ;;
        --start)        START="$2";       shift 2 ;;
        --no-moons)     NO_MOONS=true;    shift ;;
        --force)        FORCE="$2";       shift 2 ;;
        --theta)        THETA="$2";       shift 2 ;;
        --visualize)    VISUALIZE=true;   shift ;;
        --render)       RENDER=true;      shift ;;
        --fp32)         FP32=true;        shift ;;
        -h|--help)
            echo "Usage: ./run.sh [OPTIONS]"
            echo ""
            echo "Scenario selection (pick at most one):"
            echo "  --fetch              Fetch solar system IC from JPL Horizons"
            echo "  --galaxy             Generate MW-Andromeda galaxy collision IC"
            echo "  (default)            Use existing tests/sim_ic.bin"
            echo ""
            echo "Galaxy options:"
            echo "  --n N                Particles per galaxy (default: 25000)"
            echo ""
            echo "Solar options:"
            echo "  --start DATE         Fetch start date (default: 1950-01-01)"
            echo "  --no-moons           Planets only, no moons (10 bodies vs 35)"
            echo ""
            echo "Force kernel:"
            echo "  --force MODE         direct or bh (default: direct)"
            echo "  --theta T            Barnes-Hut opening angle (default: 0.5)"
            echo "  --fp32              Build using single precision (default: double)"
            echo ""
            echo "Post-processing:"
            echo "  --visualize          Open matplotlib viewer after simulation"
            echo "  --render             Open Rerun viewer after simulation"
            echo ""
            echo "  -h, --help           Show this help message"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ "$FETCH" = true ] && [ "$GALAXY" = true ]; then
    echo "ERROR: cannot combine --fetch and --galaxy"
    exit 1
fi

# Stage IC:
if [ "$FETCH" = true ]; then
    echo "Fetching solar IC from JPL Horizons..."
    MOONS_FLAG=""
    if [ "$NO_MOONS" = false ]; then MOONS_FLAG="--moons"; fi
    python3 tools/solar/fetch.py --start "$START" $MOONS_FLAG
elif [ "$GALAXY" = true ]; then
    echo "Generating MW-Andromeda galaxy IC..."
    python3 tools/galaxy/generate.py --n "$N"
fi

# Build:
echo ""
echo "Building..."
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=g++-15 \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.3/bin/nvcc \
    -DCMAKE_CUDA_HOST_COMPILER=g++-15 \
    -DFP32="$FP32" > /dev/null 2>&1
cmake --build build --target main --config Release

# Run:
echo ""
echo "Running simulation (force=$FORCE, theta=$THETA)..."
./build/main --force "$FORCE" --theta "$THETA"

# Validate (solar only):
if [ "$GALAXY" = false ]; then
    echo ""
    echo "Validating against JPL Horizons..."
    python3 tools/solar/compare.py || true
fi

# Visualize:
if [ "$VISUALIZE" = true ]; then
    echo ""
    echo "Opening matplotlib viewer..."
    if [ "$GALAXY" = true ]; then
        python3 tools/galaxy/visualize.py
    else
        python3 tools/solar/visualize.py
    fi
fi

if [ "$RENDER" = true ]; then
    echo ""
    echo "Opening Rerun viewer..."
    if [ "$GALAXY" = true ]; then
        python3 tools/galaxy/render.py
    else
        python3 tools/solar/render.py
    fi
fi

echo ""
echo "Done."
