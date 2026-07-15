@echo off
setlocal enabledelayedexpansion

echo ====================================================
echo Starting MACE-MH-1 Pre-relaxation for HAP Surfaces
echo ====================================================

REM Define paths
set MODEL_FILE=mace_pretrained_models/mace_mh_1.pt
set UTILS_DIR=../../utils

REM Ensure model is converted to .pt if needed
if not exist "%MODEL_FILE%" (
    echo Model file not found: %MODEL_FILE%
    echo Please run download_mace_model.py and convert_mace_to_mliap.py first.
    exit /b 1
)

REM Loop over all HAP_xxx_Perfect directories in 1.GeoOpt
for /d %%D in (..\HAP_*_Perfect) do (
    set DIR_NAME=%%~nxD
    echo.
    echo Processing: !DIR_NAME!
    
    set XYZ_FILE=..\!DIR_NAME!\!DIR_NAME!.xyz
    if not exist "!XYZ_FILE!" (
        REM Try Original.xyz
        set XYZ_FILE=..\!DIR_NAME!\HAP_*_Original.xyz
    )
    
    REM Conversion to LAMMPS data
    set DATA_FILE=!DIR_NAME!.data
    echo Converting XYZ to LAMMPS data...
    python %UTILS_DIR%/xyz_to_lammps.py "!XYZ_FILE!" "!DATA_FILE!"
    
    if exist "!DATA_FILE!" (
        echo Running LAMMPS Minimization via Docker...
        set DATA_FILE=!DATA_FILE!
        set OUTPUT_FILE=!DIR_NAME!_opt
        
        docker-compose run --rm -e DATA_FILE=!DATA_FILE! -e MODEL_FILE=%MODEL_FILE% -e OUTPUT_FILE=!OUTPUT_FILE! mace_prerelax
        
        if exist "!OUTPUT_FILE!.xyz" (
            echo Success for !DIR_NAME!
            move "!OUTPUT_FILE!.xyz" "..\!DIR_NAME!\!DIR_NAME!_mace_opt.xyz"
        ) else (
            echo Failed to generate output for !DIR_NAME!
        )
    ) else (
        echo Failed to create data file for !DIR_NAME!
    )
)

echo.
echo All pre-relaxations completed!
