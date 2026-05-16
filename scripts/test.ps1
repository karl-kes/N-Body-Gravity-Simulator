#Requires -Version 5.1

<#
.SYNOPSIS
    Unit Test & Benchmark Runner; builds and runs tests and optional scaling benchmark.

.EXAMPLE
    .\test.ps1                          # build and run unit tests
    .\test.ps1 -Benchmark               # also run the scaling benchmark
    .\test.ps1 -BenchOnly               # skip tests, run benchmark only
    .\test.ps1 -Benchmark -Trials 5     # pass --trials to benchmark
    .\test.ps1 -Benchmark -MaxN 8192    # pass --max-n to benchmark
#>

param(
    [switch]$Benchmark,
    [switch]$BenchOnly,
    [int]$Trials = 3,
    [int]$MaxN = 65536
)

$ErrorActionPreference = "Stop"

$RunTests = -not $BenchOnly
$RunBench = $Benchmark -or $BenchOnly

# Build
Write-Host "Building..."
cmake -B build -G "MinGW Makefiles" 2>$null | Out-Null

if ($RunTests) {
    cmake --build build --target tests --config Release
    if ($LASTEXITCODE -ne 0) { Write-Error "Build failed (tests)"; exit 1 }
}

if ($RunBench) {
    cmake --build build --target benchmark --config Release
    if ($LASTEXITCODE -ne 0) { Write-Error "Build failed (benchmark)"; exit 1 }
}

# Run tests
if ($RunTests) {
    Write-Host ""
    ./build/tests.exe
    if ($LASTEXITCODE -ne 0) { Write-Error "Tests failed"; exit 1 }
}

# Run benchmark
if ($RunBench) {
    Write-Host ""
    $BenchArgs = @("--trials", $Trials, "--max-n", $MaxN)
    ./build/benchmark.exe @BenchArgs
    if ($LASTEXITCODE -ne 0) { Write-Error "Benchmark failed"; exit 1 }
}

Write-Host ""
Write-Host "Done."
