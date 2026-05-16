#Requires -Version 5.1

<#
.SYNOPSIS
    N-Body Simulation Runner; configures, builds, and runs the full pipeline.

.EXAMPLE
    .\run.ps1                                    # defaults: dt=900, years=249, output_hours=487, with moons
    .\run.ps1 -Dt t -Years y                     # custom timestep and duration; years can be (0, 250)
    .\run.ps1 -NoMoons                           # planets only (10 bodies)
    .\run.ps1 -SkipFetch                         # reuse existing JPL data
    .\run.ps1 -Visualize                         # open Matplotlib viewer after sim
    .\run.ps1 -Render                            # open Rerun dashboard after sim
#>

param(
    [double]$Dt = 900.0,
    [int]$Years = 249,
    [int]$OutputHours = 487,
    [string]$Start = "1950-01-01",
    [switch]$NoMoons,
    [switch]$SkipFetch,
    [switch]$Visualize,
    [switch]$Render,
    [string]$Force = "direct",
    [double]$Theta = 0.5
)

$ErrorActionPreference = "Stop"
$Config = "src/Config.hpp"

# Validate:
$Remainder = [int]($OutputHours * 3600) % [int]$Dt
if ($Remainder -ne 0) {
    Write-Error "output_hours * 3600 ($OutputHours * 3600) must be divisible by dt ($Dt). Remainder: $Remainder"
    exit 1
}

# Patch Config.hpp:
Write-Host "Configuring..."
Write-Host "  dt:           $Dt s"
Write-Host "  years:        $Years"
Write-Host "  output_hours: $OutputHours"
Write-Host "  start:        $Start"
Write-Host "  moons:        $(if ($NoMoons) { 'no' } else { 'yes' })"

$Content = Get-Content $Config -Raw
$Content = $Content -replace '(inline static constexpr double dt\{)[^}]*', "`${1} $Dt "
$Content = $Content -replace '(inline static constexpr std::size_t num_years\{)[^}]*', "`${1} $Years "
$Content = $Content -replace '(inline static constexpr std::size_t output_hours\{)[^}]*', "`${1} $OutputHours "
Set-Content $Config -Value $Content -NoNewline

Write-Host "  Config.hpp updated."

# Fetch JPL data:
if (-not $SkipFetch) {
    Write-Host ""
    Write-Host "Fetching JPL Horizons data..."
    $FetchArgs = @("src/jpl_compare.py", "fetch", "--start", $Start)
    if (-not $NoMoons) { $FetchArgs += "--moons" }
    python $FetchArgs
    if ($LASTEXITCODE -ne 0) { Write-Error "JPL fetch failed"; exit 1 }
} else {
    Write-Host ""
    Write-Host "Skipping JPL fetch (-SkipFetch)"
}

# Build:
Write-Host ""
Write-Host "Building..."

cmake -B build -G "MinGW Makefiles" 2>$null | Out-Null
cmake --build build --config Release
if ($LASTEXITCODE -ne 0) { Write-Error "Build failed"; exit 1 }

# Run:
Write-Host ""
Write-Host "Running simulation..."

./build/main.exe --force $Force --theta $Theta
if ($LASTEXITCODE -ne 0) { Write-Error "Simulation failed"; exit 1 }

# Validate:
Write-Host ""
Write-Host "Validating against JPL Horizons..."

python src/jpl_compare.py compare
if ($LASTEXITCODE -ne 0) { Write-Error "Validation failed"; exit 1 }

# Visualize:
if ($Visualize) {
    Write-Host ""
    python src/visualize.py
}

if ($Render) {
    Write-Host ""
    python src/render.py
}

Write-Host ""
Write-Host "Done."
