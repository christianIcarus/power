---
name: run-solar-efficiency
description: Build, run, and drive solar_efficiency.py (the PX4 flight-log solar-array efficiency analyzer). Use when asked to run this tool, test it, verify it still works, or generate/screenshot its efficiency plots from a .ulg log.
---

This repo is a single Python CLI script, [solar_efficiency.py](../../../solar_efficiency.py), that reads a PX4 `.ulg` flight log and estimates pre-MPPT solar-array efficiency against a modeled clear-sky irradiance. Drive it via the smoke driver at `.claude/skills/run-solar-efficiency/driver.ps1` (Windows PowerShell) - it resolves Python, installs deps if missing, runs the script against a real log, and verifies the CSV/plots it produces are non-trivial, not just that the process exited 0.

There is no fixture log committed to the repo - real `.ulg` logs are large (the one used to verify this skill is 1.1 GB) and machine-specific. You must supply a real log path.

## Prerequisites

Windows only - the script's paths and this driver assume Windows (e.g. `open_in_vscode` shells out to `code`, the documented default log path is a `C:\Users\...` path). No OS packages needed beyond Python itself.

```powershell
winget install -e --id Python.Python.3.12 --accept-package-agreements --accept-source-agreements --silent
```

This installs to `%LOCALAPPDATA%\Programs\Python\Python312\python.exe`. **Do not rely on `python`/`py` on `PATH`** - see Gotchas. The driver script resolves the real interpreter itself; you don't need to do this manually.

## Setup

```powershell
$py = "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe"
& $py -m pip install -r requirements.txt
```

`requirements.txt` covers `pyulog`, `pvlib`, `pandas`, `numpy`, `matplotlib`, `timezonefinder`. (The driver does this automatically if imports fail - manual install is only needed if you're not using the driver.)

## Build

No build step - it's a plain script, run directly.

## Run (agent path)

Use the driver. Fast smoke mode (`--no-plot`, skips PNG rendering, ~30s against a 1.1 GB / 2h12m log) is the default - it still parses the full log, runs the irradiance/efficiency model, and writes the CSV, which is enough to prove the pipeline works end to end:

```powershell
.\.claude\skills\run-solar-efficiency\driver.ps1 -Ulog "C:\path\to\flight.ulg"
```

Pass `-Full` to also render all 6 plots (main summary, per-string, POA-per-string, string-1, string-1 %-diff, panel-normal-angle sweep) - takes longer (~55s on the 1.1 GB verification log, vs. ~30s for the fast path):

```powershell
.\.claude\skills\run-solar-efficiency\driver.ps1 -Ulog "C:\path\to\flight.ulg" -Full
```

Optionally pass `-OutputDir <path>` (defaults to a `smoke_out` folder next to the log). The driver prints `SMOKE OK` and exits 0 only if the CSV (and, under `-Full`, the main PNG) actually exist and are a plausible size - not just that the Python process returned 0.

To *look at* a plot after a `-Full` run (e.g. to visually confirm it's not a blank frame), read the PNG with an image-capable tool - it's a static file, no display needed:

```
<OutputDir>\<ulog-stem>_solar_efficiency.png
```

## Run (human path)

```powershell
python solar_efficiency.py --ulog "C:\path\to\flight.ulg"
```

Without `--no-open`, each generated plot is shelled out to `code <file>.png` (VS Code) as it's produced. Without `--output-dir`, output lands next to the log file. Ctrl-C is not needed - the script runs to completion and exits.

## Test

There is no test suite in this repo (no `pytest`/`unittest` files) - the driver's smoke run against a real log is the verification path.

---

## Gotchas

- **`python`/`py` on `PATH` may be Windows Store stub aliases, not real Python.** On a machine where Python was never installed, `python.exe`/`python3.exe` under `AppData\Local\Microsoft\WindowsApps` exist as app-execution-alias stubs that just open the Microsoft Store when run - they are not an interpreter. `py.exe` may not exist at all. The driver resolves this by checking `Get-Command python` and rejecting anything under `WindowsApps`, falling back to the known `Programs\Python\Python312\python.exe` install path.
- **PATH changes don't persist across separate shell invocations in this harness.** Each tool call is a fresh process; `winget install`'s PATH update isn't visible until a new shell/session. The driver sidesteps this by resolving Python's absolute path itself rather than trusting `PATH`.
- **`timezonefinder` is a hard runtime dependency missing from `requirements.txt`** (imported at module load, `from timezonefinder import TimezoneFinder`) - already fixed in this repo's `requirements.txt`, but if you see `ModuleNotFoundError: No module named 'timezonefinder'` on an older checkout, `pip install timezonefinder` directly.
- **The script's own default `--ulog` path is hardcoded to a specific machine/date** (`solar_efficiency.py`'s `DEFAULT_ULOG`) and will not exist elsewhere - always pass `--ulog` explicitly.
- **Don't pipe the script's stderr through `2>&1` in Windows PowerShell 5.1.** matplotlib's one-time "building the font cache" message (and any other stderr chatter) gets wrapped in a `NativeCommandError` and flips `$?` to `$false` even when the process exit code is 0. Check `$LASTEXITCODE`, not `$?` - the driver does this.
- **A `pv_power_w_0 contributed only 0.00%` warning is a data-content signal, not a driver failure** - it means one MPPT string reads ~0 for that particular flight (e.g. disconnected string), and the script says so itself; it still exits 0 and produces valid output for the string that is connected.

## Troubleshooting

- **`ModuleNotFoundError: No module named 'timezonefinder'`**: not in `requirements.txt` on older checkouts (fixed now). `pip install timezonefinder`.
- **`python : The term 'python' is not recognized...` (PowerShell) even right after installing Python**: the install succeeded but the current shell's `PATH` is stale. Use the absolute path (`$env:LOCALAPPDATA\Programs\Python\Python312\python.exe`) instead of relying on `PATH`, as the driver does.
- **`Python was not found; run without arguments to install from the Microsoft Store...`**: you invoked the Windows Store alias stub, not real Python. Install via winget (see Prerequisites) and use the absolute path to the real `python.exe`.
