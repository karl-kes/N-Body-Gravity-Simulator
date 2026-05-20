#Requires -Version 5.1

<#
.SYNOPSIS
    Simulation runner. Builds and runs main against tests/sim_ic.bin,
    which is produced by jpl_compare.py fetch (solar) or
    generate_galaxies.py (galaxy).

.EXAMPLE
    .\scripts\run.ps1                              # use existing sim_ic.bin
    .\scripts\run.ps1 -Fetch                       # fetch solar IC, then run
    .\scripts\run.ps1 -Galaxy -N 25000             # generate galaxy IC, then run
    .\scripts\run.ps1 -Force bh -Theta 0.5         # pass through to main
    .\scripts\run.ps1 -Visualize                   # matplotlib viewer after sim
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
    Write-Error "cannot combine -Fetch and -Galaxy"
    exit 1
}

if ($Fetch) {
    Write-Host "Fetching solar IC..."
    $FetchArgs = @("src/jpl_compare.py", "fetch", "--start", $Start)
    if (-not $NoMoons) { $FetchArgs += "--moons" }
    python $FetchArgs
    if ($LASTEXITCODE -ne 0) { Write-Error "JPL fetch failed"; exit 1 }
} elseif ($Galaxy) {
    Write-Host "Generating galaxy IC..."
    python src/generate_galaxies.py --n $N
    if ($LASTEXITCODE -ne 0) { Write-Error "Galaxy generation failed"; exit 1 }
}

Write-Host ""
Write-Host "Building..."
cmake -B build -G "MinGW Makefiles" 2>$null | Out-Null
cmake --build build --config Release
if ($LASTEXITCODE -ne 0) { Write-Error "Build failed"; exit 1 }

Write-Host ""
Write-Host "Running simulation (force=$Force, theta=$Theta)..."
./build/main.exe --force $Force --theta $Theta
if ($LASTEXITCODE -ne 0) { Write-Error "Simulation failed"; exit 1 }

if (-not $Galaxy) {
    Write-Host ""
    Write-Host "Validating against JPL Horizons..."
    python src/jpl_compare.py compare
}

if ($Visualize) {
    Write-Host ""
    if ($Galaxy) {
        python src/visualize_galaxies.py
    } else {
        python src/visualize.py
    }
}

if ($Render) {
    Write-Host ""
    python src/render.py
}

Write-Host ""
Write-Host "Done."