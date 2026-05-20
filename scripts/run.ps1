#Requires -Version 5.1

<#
.SYNOPSIS
    Simulation Pipeline Runner; generates the IC, builds, runs, validates, and visualizes.

.EXAMPLE
    .\run.ps1                                  # use existing tests/sim_ic.bin
    .\run.ps1 -Fetch                           # fetch solar IC from JPL, then run
    .\run.ps1 -Galaxy                          # generate MW-Andromeda IC, then run
    .\run.ps1 -Galaxy -N 50000                 # 50k particles per galaxy
    .\run.ps1 -Force bh -Theta 0.5             # Barnes-Hut force at theta = 0.5
    .\run.ps1 -Galaxy -Force bh -Visualize
    .\run.ps1 -Fetch -Visualize -Render        # full solar pipeline with both viewers
#>

param(
    [switch]$Fetch,
    [switch]$Galaxy,
    [int]$N = 25000,
    [string]$Start = "1950-01-01",
    [switch]$NoMoons,
    [string]$Force = "direct",
    [double]$Theta = 0.5,
    [switch]$Visualize,
    [switch]$Render
)

$ErrorActionPreference = "Stop"

if ($Fetch -and $Galaxy) {
    Write-Error "ERROR: cannot combine -Fetch and -Galaxy"
    exit 1
}

# Stage IC:
if ($Fetch) {
    Write-Host "Fetching solar IC from JPL Horizons..."
    $FetchArgs = @("tools/solar/fetch.py", "--start", $Start)
    if (-not $NoMoons) { $FetchArgs += "--moons" }
    python $FetchArgs
    if ($LASTEXITCODE -ne 0) { Write-Error "JPL fetch failed"; exit 1 }
} elseif ($Galaxy) {
    Write-Host "Generating MW-Andromeda galaxy IC..."
    python tools/galaxy/generate.py --n $N
    if ($LASTEXITCODE -ne 0) { Write-Error "Galaxy generation failed"; exit 1 }
}

# Build:
Write-Host ""
Write-Host "Building..."
cmake -B build -G "MinGW Makefiles" 2>$null | Out-Null
cmake --build build --target main --config Release
if ($LASTEXITCODE -ne 0) { Write-Error "Build failed"; exit 1 }

# Run:
Write-Host ""
Write-Host "Running simulation (force=$Force, theta=$Theta)..."
./build/main.exe --force $Force --theta $Theta
if ($LASTEXITCODE -ne 0) { Write-Error "Simulation failed"; exit 1 }

# Validate (solar only):
if (-not $Galaxy) {
    Write-Host ""
    Write-Host "Validating against JPL Horizons..."
    python tools/solar/compare.py
}

# Visualize:
if ($Visualize) {
    Write-Host ""
    Write-Host "Opening matplotlib viewer..."
    if ($Galaxy) {
        python tools/galaxy/visualize.py
    } else {
        python tools/solar/visualize.py
    }
}

if ($Render) {
    Write-Host ""
    Write-Host "Opening Rerun viewer..."
    if ($Galaxy) {
        python tools/galaxy/render.py
    } else {
        python tools/solar/render.py
    }
}

Write-Host ""
Write-Host "Done."