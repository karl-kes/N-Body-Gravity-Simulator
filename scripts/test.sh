#!/usr/bin/env bash
set -euo pipefail

# Unit Test & Benchmark Runner
# Builds and runs the test suite and optional scaling benchmark.
#
# Usage:
#   ./test.sh                               # build and run unit tests
#   ./test.sh --benchmark                   # also run the scaling benchmark
#   ./test.sh --bench-only                  # skip tests, run benchmark only
#   ./test.sh --benchmark --trials 5        # pass --trials to benchmark
#   ./test.sh --benchmark --max-n 8192      # pass --max-n to benchmark

BENCHMARK=false
BENCH_ONLY=false
TRIALS=3
MAX_N=65536

while [[ $# -gt 0 ]]; do
    case $1 in
        --benchmark)    BENCHMARK=true;   shift ;;
        --bench-only)   BENCH_ONLY=true;  shift ;;
        --trials)       TRIALS="$2";      shift 2 ;;
        --max-n)        MAX_N="$2";       shift 2 ;;
        -h|--help)
            echo "Usage: ./test.sh [OPTIONS]"
            echo ""
            echo "Mode:"
            echo "  (default)            Build and run unit tests"
            echo "  --benchmark          Also run the scaling benchmark after tests"
            echo "  --bench-only         Skip tests, run benchmark only"
            echo ""
            echo "Benchmark options:"
            echo "  --trials N           Trials per N value (default: 3)"
            echo "  --max-n N            Maximum N for benchmark sweep (default: 65536)"
            echo ""
            echo "  -h, --help           Show this help message"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

RUN_TESTS=true
RUN_BENCH=false
if [ "$BENCH_ONLY" = true ]; then RUN_TESTS=false; RUN_BENCH=true; fi
if [ "$BENCHMARK" = true ];  then RUN_BENCH=true; fi

# Build:
echo "Building..."
cmake -B build -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1

if [ "$RUN_TESTS" = true ]; then
    cmake --build build --target tests --config Release
fi

if [ "$RUN_BENCH" = true ]; then
    cmake --build build --target benchmark --config Release
fi

# Run tests:
if [ "$RUN_TESTS" = true ]; then
    echo ""
    echo "Running unit tests..."
    ./build/tests
fi

# Run benchmark:
if [ "$RUN_BENCH" = true ]; then
    echo ""
    echo "Running scaling benchmark..."
    ./build/benchmark --trials "$TRIALS" --max-n "$MAX_N"
fi

echo ""
echo "Done."