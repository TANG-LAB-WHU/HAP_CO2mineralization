@echo off
REM ====================================================================
REM  Windows batch script to generate a clean XYZ, then launch a
REM  GPU-accelerated CP2K geometry optimisation for HAP002_THICKNESS1
REM ====================================================================

REM --------------------------------------------------------------------
REM 0. Detect hkl facet from folder name and set PROJECT_BASE
REM    Expected folder names like: HAP_112_Perfect, HAP_211_*, etc.
REM    Fallback to HAP_hkl if not detected
REM --------------------------------------------------------------------
setlocal EnableDelayedExpansion
set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR:~0,-1%") do set "FOLDER_NAME=%%~nxI"

set "PROJECT_BASE=HAP_hkl"
for /f "tokens=1,2 delims=_" %%A in ("!FOLDER_NAME!") do (
    if /I "%%A"=="HAP" (
        if not "%%B"=="" set "PROJECT_BASE=HAP_%%B"
    )
)
echo Detected PROJECT_BASE: %PROJECT_BASE%

REM --------------------------------------------------------------------
REM 1. Generate clean structure (requires Python + ASE installed)
REM --------------------------------------------------------------------
echo Generating clean structure with Python...
set "PROJECT_BASE=%PROJECT_BASE%"
python "%~dp0remove_duplicate_atoms.py"
if %errorlevel% neq 0 (
    echo ERROR: remove_duplicate_atoms.py failed. Aborting.
    pause
    exit /b 1
)

REM --------------------------------------------------------------------
REM 2. Prepare run directory and facet-specific CP2K input
REM --------------------------------------------------------------------
if not exist "cp2k_run" mkdir cp2k_run
copy /Y "%PROJECT_BASE%_clean.xyz" "cp2k_run\%PROJECT_BASE%_clean.xyz" > nul
copy /Y "cp2k_run\HAP_hkl_opt.inp" "cp2k_run\%PROJECT_BASE%_opt.inp" > nul
powershell -NoProfile -ExecutionPolicy Bypass -Command "(Get-Content 'cp2k_run\%PROJECT_BASE%_opt.inp') -replace '@SET PROJECT HAP_hkl','@SET PROJECT %PROJECT_BASE%' | Set-Content -Encoding ASCII 'cp2k_run\%PROJECT_BASE%_opt.inp'"

REM --------------------------------------------------------------------
REM 3. Launch CP2K in Docker with facet-specific filenames
REM    docker compose will substitute ${PROJECT_BASE} from this environment
REM --------------------------------------------------------------------
docker compose -f docker-compose-cp2k.yml down --remove-orphans
docker compose -f docker-compose-cp2k.yml up --abort-on-container-exit

pause 