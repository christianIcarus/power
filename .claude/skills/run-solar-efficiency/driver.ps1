#Requires -Version 5.1
<#
.SYNOPSIS
  Smoke-drives solar_efficiency.py: resolves a working Python, installs
  deps if missing, runs the analysis against a real .ulg log, and checks
  that it actually produced output (not just exit-code 0).

.PARAMETER Ulog
  Path to a PX4 .ulg flight log. Required - there is no fixture log
  committed to the repo (real logs are large and machine-specific).

.PARAMETER OutputDir
  Where the CSV/PNGs land. Defaults to a "smoke_out" folder next to the
  log file.

.PARAMETER Full
  Also generate all plots (default is --no-plot: parses the log, runs
  the irradiance/efficiency model, writes the CSV, skips PNG rendering -
  much faster and enough to prove the pipeline works end to end).

.EXAMPLE
  .\.claude\skills\run-solar-efficiency\driver.ps1 -Ulog "C:\Users\me\Downloads\flight.ulg"

.EXAMPLE
  .\.claude\skills\run-solar-efficiency\driver.ps1 -Ulog "C:\Users\me\Downloads\flight.ulg" -Full
#>
param(
    [Parameter(Mandatory = $true)][string]$Ulog,
    [string]$OutputDir,
    [switch]$Full
)

$ErrorActionPreference = "Stop"
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..\..")
$script = Join-Path $repoRoot "solar_efficiency.py"

if (-not (Test-Path $Ulog)) {
    Write-Error "Ulog not found: $Ulog"
    exit 1
}
if (-not $OutputDir) {
    $OutputDir = Join-Path (Split-Path $Ulog -Parent) "smoke_out"
}

# --- Resolve a real Python (the Windows Store 'python'/'py' aliases are stubs
#     that open the Store instead of running anything - see Gotchas). ---
function Resolve-Python {
    $candidates = @()
    $onPath = Get-Command python -ErrorAction SilentlyContinue
    if ($onPath -and $onPath.Source -notmatch "WindowsApps") { $candidates += $onPath.Source }
    $candidates += "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe"
    $candidates += "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe"
    foreach ($c in $candidates) {
        if ($c -and (Test-Path $c)) { return $c }
    }
    return $null
}

$py = Resolve-Python
if (-not $py) {
    Write-Host "No usable Python found - installing Python 3.12 via winget ..."
    winget install -e --id Python.Python.3.12 --accept-package-agreements --accept-source-agreements --silent
    $py = Resolve-Python
    if (-not $py) { Write-Error "Python install did not produce a usable python.exe"; exit 1 }
}
Write-Host "Using Python: $py ($(& $py --version))"

# --- Verify deps import; install requirements.txt (+ timezonefinder) if not. ---
& $py -c "import pvlib, pandas, numpy, matplotlib, pyulog, timezonefinder" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Missing dependencies - installing from requirements.txt ..."
    & $py -m pip install -r (Join-Path $repoRoot "requirements.txt")
    if ($LASTEXITCODE -ne 0) { Write-Error "pip install -r requirements.txt failed"; exit 1 }
    & $py -c "import pvlib, pandas, numpy, matplotlib, pyulog, timezonefinder"
    if ($LASTEXITCODE -ne 0) { Write-Error "Dependencies still missing after install"; exit 1 }
}

# --- Run the analysis. Capture stdout/stderr SEPARATELY - piping a native
#     exe's stderr through 2>&1 in Windows PowerShell 5.1 wraps each line in
#     a NativeCommandError and flips $? to false even on exit code 0. Check
#     $LASTEXITCODE, not $?. ---
$argList = @("--ulog", $Ulog, "--output-dir", $OutputDir, "--no-open")
if (-not $Full) { $argList += "--no-plot" }

Write-Host "Running: python solar_efficiency.py $($argList -join ' ')"
$stdout = & $py $script @argList
$exitCode = $LASTEXITCODE
$stdout | ForEach-Object { Write-Host $_ }

if ($exitCode -ne 0) {
    Write-Error "solar_efficiency.py exited $exitCode"
    exit $exitCode
}

# --- Prove it actually produced output, not just a clean exit. ---
$stem = [System.IO.Path]::GetFileNameWithoutExtension($Ulog)
$csv = Join-Path $OutputDir "${stem}_solar_efficiency.csv"
$ok = $true

if (-not (Test-Path $csv) -or (Get-Item $csv).Length -lt 1000) {
    Write-Host "FAIL: CSV missing or suspiciously small: $csv"
    $ok = $false
} else {
    $rows = (Get-Content $csv | Measure-Object -Line).Lines
    Write-Host "OK: CSV written ($rows lines) -> $csv"
}

if ($Full) {
    $plot = Join-Path $OutputDir "${stem}_solar_efficiency.png"
    if (-not (Test-Path $plot) -or (Get-Item $plot).Length -lt 10000) {
        Write-Host "FAIL: main plot missing or suspiciously small: $plot"
        $ok = $false
    } else {
        Write-Host "OK: plot written ($((Get-Item $plot).Length) bytes) -> $plot"
    }
}

if (-not $ok) { exit 1 }
Write-Host "SMOKE OK"
exit 0
