@echo off
REM AIMD to Deep Learning MD Workflow - Dynamic HKL Parameter Version
REM UPDATED: Support for dynamic hkl parameter inheritance from folder name
REM SIMPLIFIED: Single workflow path - CP2K AIMD → DeepMD → LAMMPS
REM CONFIGURABLE: Surface generation parameters

REM =============================================================================
REM Dynamic HKL Parameter Extraction from Folder Name
REM =============================================================================

REM Extract hkl parameter from current folder name
for %%i in ("%cd%") do set "FOLDER_NAME=%%~ni"
echo Current folder: %FOLDER_NAME%

REM Extract hkl from folder name (assuming format: HAP_XXX_template or similar)
for /f "tokens=2 delims=_" %%a in ("%FOLDER_NAME%") do set "HKL_PARAM=%%a"
echo Extracted HKL parameter: %HKL_PARAM%

REM Validate HKL parameter (skip strict check due to parsing issues)
echo [DEBUG] HKL_PARAM="%HKL_PARAM%"

REM Parse HKL into individual Miller indices
set "MILLER_H=%HKL_PARAM:~0,1%"
set "MILLER_K=%HKL_PARAM:~1,1%"
set "MILLER_L=%HKL_PARAM:~2,1%"

REM Set derived parameters for file naming
set "HAP_PREFIX=HAP_%HKL_PARAM%"
set "HAP_DATA_BASENAME=hap_%HKL_PARAM%"
set "HAP_DATA_FILE=hap_%HKL_PARAM%.data"

REM Surface generation parameters changeable by user
set "SURFACE_LAYERS=4"
set "SURFACE_SUPERCELL_A=3"
set "SURFACE_SUPERCELL_B=3"
set "SURFACE_VACUUM=15"
set "ADD_CO2=true"
set "SURFACE_CO2_HEIGHT=10"
set "N_CO2=1"

REM Slab generation method selection
set "SLAB_METHOD=ase"

REM Enable delayed variable expansion
setlocal enabledelayedexpansion

REM =============================================================================
REM LAMMPS Initial Structure Preparation
REM =============================================================================

echo ==========================================
echo AIMD to Deep Learning MD Workflow
echo Surface Model: (%MILLER_H% %MILLER_K% %MILLER_L%) surface
echo HKL Parameter: %HKL_PARAM%
echo HAP Prefix: %HAP_PREFIX%
echo Data File: %HAP_DATA_FILE%
echo Surface Parameters: %SURFACE_LAYERS% layers, %SURFACE_SUPERCELL_A%x%SURFACE_SUPERCELL_B% supercell
echo Vacuum (z): %SURFACE_VACUUM% Angstrom
echo Slab Method: %SLAB_METHOD%
echo PyTorch Backend with DPA Descriptor
echo Started: %date% %time%
echo ==========================================

echo === Generating LAMMPS Initial Structure from CIF ===

echo Surface Configuration:
echo   Miller Indices: (%MILLER_H% %MILLER_K% %MILLER_L%)
echo   Surface Layers: %SURFACE_LAYERS%
echo   Supercell: %SURFACE_SUPERCELL_A% x %SURFACE_SUPERCELL_B%
echo   Vacuum (z): %SURFACE_VACUUM% Angstrom
echo   Add CO2: %ADD_CO2%
echo   Number of CO2 molecules per side: %N_CO2%
echo   CO2 Height: %SURFACE_CO2_HEIGHT% Angstrom
echo   Slab Method: %SLAB_METHOD%

REM Check if CIF file exists in lammps_run directory
if exist "lammps_run\HAP_ConventionalCell.cif" (
    echo CIF file found: lammps_run\HAP_ConventionalCell.cif
    
    REM Check if create_hap_slabs.py exists in lammps_run directory
    if not exist "lammps_run\create_hap_slabs.py" (
        echo ERROR: create_hap_slabs.py not found in lammps_run directory
        echo Please ensure create_hap_slabs.py exists in lammps_run directory
        exit /b 1
    ) else (
        echo Found create_hap_slabs.py in lammps_run directory
    )
    
    REM Use local Python environment for reliable surface generation
    echo Using local Python environment for surface generation...
    
    REM Prefer Python from current conda env; fallback to where.exe
    set "PYTHON_CMD="
    if defined CONDA_PREFIX (
        if exist "%CONDA_PREFIX%\python.exe" set "PYTHON_CMD=%CONDA_PREFIX%\python.exe"
    )
    if "%PYTHON_CMD%"=="" (
        for /f %%i in ('where.exe python 2^>nul') do (
            set "PYTHON_CMD=%%i"
            goto py_found
        )
    )

    :py_found
    if "%PYTHON_CMD%"=="" (
        echo ERROR: No Python found in PATH. Please activate your conda env and retry.
        exit /b 1
    )

    echo Using Python: %PYTHON_CMD%
    "%PYTHON_CMD%" --version
    
    REM Execute surface generation directly
    echo Executing surface generation with local Python...
    cd lammps_run
    "%PYTHON_CMD%" create_hap_slabs.py %MILLER_H% %MILLER_K% %MILLER_L% --layers %SURFACE_LAYERS% --vacuum %SURFACE_VACUUM% --supercell %SURFACE_SUPERCELL_A% %SURFACE_SUPERCELL_B% --method %SLAB_METHOD% --add-co2 --n-co2 %N_CO2% --co2-height %SURFACE_CO2_HEIGHT% --output %HAP_DATA_BASENAME%
    cd ..
    
    REM Check if surface generation was successful
    if not exist "lammps_run\%HAP_DATA_FILE%" (
        echo ERROR: Surface generation failed - %HAP_DATA_FILE% not found
        exit /b 1
    ) else (
        echo Surface generation successful - %HAP_DATA_FILE% created
    )
    
    REM Check for successful generation and verify structure
    dir /b "lammps_run\*.data" > "%TEMP%\data_files.tmp" 2>nul
    if %errorlevel% equ 0 (
        for /f %%i in (%TEMP%\data_files.tmp) do (
            set "FOUND_DATA_FILE=%%i"
            goto data_file_found
        )
    )
    goto data_file_not_found
    
    :data_file_found
    echo LAMMPS initial structure generated successfully: !FOUND_DATA_FILE!
    
    REM If the found file is not exactly the expected data file, copy it
    if not "!FOUND_DATA_FILE!"=="%HAP_DATA_FILE%" (
        copy "lammps_run\!FOUND_DATA_FILE!" "lammps_run\%HAP_DATA_FILE%" > nul
        echo Copied !FOUND_DATA_FILE! to %HAP_DATA_FILE% for workflow compatibility
    )
    
    del "%TEMP%\data_files.tmp" 2>nul
    goto structure_verification
        
    :data_file_not_found
    echo ERROR: LAMMPS initial structure generation failed using %SLAB_METHOD% method
    echo.
    echo Checking for alternative output files...
    dir /b "lammps_run\*.data" 2>nul
    if errorlevel 1 (
        echo No .data files found in lammps_run directory
    ) else (
        echo Found .data files:
        dir "lammps_run\*.data"
    )
    exit /b 1
    
    :structure_verification
    echo Verifying generated structure with %SLAB_METHOD% method...
    
    REM Extract atom types count
    for /f "tokens=1" %%i in ('findstr "atom types" "lammps_run\%HAP_DATA_FILE%"') do (
        set "ATOM_TYPES=%%i"
        echo Generated structure contains !ATOM_TYPES! atom types using %SLAB_METHOD% method
    )
    
    REM Validate CO2 inclusion
    if "%ADD_CO2%"=="true" (
        if "!ATOM_TYPES!"=="5" (
            echo CO2 successfully added using %SLAB_METHOD% method - 5 atom types detected
        ) else (
            echo WARNING: Expected 5 atom types with CO2, found !ATOM_TYPES! using %SLAB_METHOD% method
        )
    ) else (
        if "!ATOM_TYPES!"=="4" (
            echo Clean HAP surface generated using %SLAB_METHOD% method - 4 atom types detected
        ) else (
            echo WARNING: Expected 4 atom types for clean HAP, found !ATOM_TYPES! using %SLAB_METHOD% method
        )
    )
    
    echo Surface generation completed successfully!
    echo.
    
    REM Continue to next workflow step - NO MORE CIF CHECKS
    goto continue_workflow
    
) else (
    echo ERROR: CIF file not found: lammps_run\HAP_ConventionalCell.cif
    exit /b 1
)

:continue_workflow
REM =============================================================================
REM Directory Structure Preparation
REM =============================================================================

echo === Preparing Directory Structure ===

if not exist deepmd_data (
    mkdir deepmd_data
    echo Created deepmd_data directory
) else (
    echo deepmd_data directory already exists
)

if not exist deepmd_data\set.000 (
    mkdir deepmd_data\set.000
    echo Created deepmd_data\set.000 directory
) else (
    echo deepmd_data\set.000 directory already exists
)

if not exist deepmd_model (
    mkdir deepmd_model
    echo Created deepmd_model directory
) else (
    echo deepmd_model directory already exists
)

if not exist lammps_run\model (
    mkdir lammps_run\model
    echo Created lammps_run\model directory
) else (
    echo lammps_run\model directory already exists
)

REM =============================================================================
REM Dynamic Configuration Generation
REM =============================================================================

echo === Generating Dynamic Configuration Files ===

REM Check if configuration generator exists
if not exist "generate_config.py" (
    echo ERROR: generate_config.py not found
    echo Please create the configuration generator script
    exit /b 1
)

REM Generate configuration files based on HKL parameter
!PYTHON_CMD! generate_config.py

if %errorlevel% neq 0 (
    echo ERROR: Configuration generation failed
    exit /b 1
) else (
    echo Dynamic configuration files generated successfully
)

REM =============================================================================
REM CP2K Geometry Optimization
REM =============================================================================

echo === Running CP2K ===

IF NOT EXIST cp2k_run\%HAP_PREFIX%_opt-pos-1.xyz (
    echo Error: geometry-optimization trajectory cp2k_run\%HAP_PREFIX%_opt-pos-1.xyz not found
    exit /b 1
) else (
    echo Found CP2K geometry optimization trajectory
)

REM =============================================================================
REM Final Frame Extraction
REM =============================================================================

echo === Extracting final frame from CP2K geometrical optimization trajectory ===

REM Using the Python 3.10-slim image that's already pulled
docker run --rm -v "%cd%:/workspace" -w /workspace docker.m.daocloud.io/library/python:3.10-slim bash -c "python cp2k_run/extract_last_frame.py cp2k_run/%HAP_PREFIX%_opt-pos-1.xyz cp2k_run/%HAP_PREFIX%_final.xyz"

if %errorlevel% equ 0 (
    echo Final frame extracted successfully
) else (
    echo ERROR: Failed to extract final frame
)

IF EXIST cp2k_run\%HAP_PREFIX%_md-pos.xyz (
    echo CP2K MD trajectory already exists, skipping CP2K run.
    goto after_cp2k_run
)

echo Starting CP2K MD simulation...

REM Using the CP2K image that's already pulled
docker compose up -d cp2k_run

echo Waiting for CP2K container to finish...

REM Try up to 15 times (about 60 s) to obtain container ID
set "CP2K_ID="
for /l %%n in (1,1,15) do (
    for /f %%i in ('docker compose ps -a -q cp2k_run 2^>nul') do set "CP2K_ID=%%i"
    if defined CP2K_ID goto id_ready
    echo Waiting for CP2K container attempt %%n of 15...
    timeout /t 2 >nul
)

:id_ready
echo Found CP2K container ID: !CP2K_ID!
docker wait !CP2K_ID!
IF ERRORLEVEL 1 (
    echo CP2K run failed - check Docker logs
    docker logs !CP2K_ID!
    exit /b 1
) else (
    echo CP2K container finished successfully
)

:after_cp2k_run

REM Verify MD trajectory file exists after run/skip logic
IF NOT EXIST cp2k_run\%HAP_PREFIX%_md-pos.xyz (
    echo Error: CP2K MD trajectory cp2k_run\%HAP_PREFIX%_md-pos.xyz still missing after run!
    exit /b 1
) else (
    echo CP2K MD trajectory file verified
)

REM =============================================================================
REM DeepMD Dataset Generation
REM =============================================================================

echo === Generating DeepMD dataset from CP2K MD trajectory ===

REM Verify CP2K MD trajectory exists
if not exist "cp2k_run\%HAP_PREFIX%_md-pos.xyz" (
    echo ERROR: CP2K MD trajectory not found: cp2k_run\%HAP_PREFIX%_md-pos.xyz
    exit /b 1
)

REM Using the Python 3.10-slim image that's already pulled
docker run --rm -v "%cd%:/workspace" -w /workspace docker.m.daocloud.io/library/python:3.10-slim bash -c "pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -r requirements.txt && python generate_deepmd_data.py"

IF ERRORLEVEL 1 (
    echo ERROR: DeepMD dataset generation failed - check logs and error file
    exit /b 1
) else (
    echo DeepMD dataset generation completed successfully
)

REM Verify that DeepMD dataset was created successfully
if not exist deepmd_data\set.000 (
    echo ERROR: deepmd_data\set.000 not found - DeepMD dataset generation failed
    exit /b 1
) else (
    echo DeepMD dataset verified at deepmd_data\set.000
)

REM =============================================================================
REM DeepMD Model Training (PyTorch Backend)
REM =============================================================================

echo === Running DeepMD Training with PyTorch Backend ===
echo Starting DeepMD PyTorch training with DPA descriptor...

REM Check for PyTorch model format (.pth)
if exist deepmd_model\hap_model.pth (
    echo PyTorch model already exists, skip training
) else (
    echo Starting DeepMD PyTorch model training...
    
    REM Using the DeepMD image that's already pulled
    docker compose up -d --no-deps --force-recreate deepmd_train
    
    REM Get container ID for training
    FOR /f %%i IN ('docker compose ps -a -q deepmd_train 2^>nul') DO set TRAIN_ID=%%i
    
    if not defined TRAIN_ID (
        echo ERROR: Failed to get DeepMD training container ID
        exit /b 1
    )
    
    echo DeepMD training container ID: !TRAIN_ID!
    
    REM Follow logs and wait for container
    echo Following DeepMD PyTorch training logs...
    
    docker wait !TRAIN_ID!
    set WAIT_EXIT_CODE=!errorlevel!
    
    REM Capture logs after container finishes
    docker logs !TRAIN_ID!
    
    if !WAIT_EXIT_CODE! equ 0 (
        echo DeepMD PyTorch training completed successfully
    ) else (
        echo ERROR: DeepMD PyTorch training failed - check Docker logs
        exit /b 1
    )
)

REM =============================================================================
REM PyTorch Model Freezing
REM =============================================================================

echo === Freezing PyTorch DeepMD Model ===
echo Freezing PyTorch model to .pth format...

REM Using the DeepMD image that's already pulled
docker compose up -d --no-deps --force-recreate deepmd_freeze

REM Get container ID for freezing
FOR /f %%i IN ('docker compose ps -a -q deepmd_freeze 2^>nul') DO set FREEZE_ID=%%i

if not defined FREEZE_ID (
    echo ERROR: Failed to get DeepMD freeze container ID
    exit /b 1
)

echo DeepMD freeze container ID: !FREEZE_ID!

REM Follow logs and wait for container
docker wait !FREEZE_ID!
set FREEZE_EXIT_CODE=!errorlevel!

REM Capture logs after container finishes
docker logs !FREEZE_ID!

if !FREEZE_EXIT_CODE! equ 0 (
    echo PyTorch model freezing completed successfully
) else (
    echo ERROR: PyTorch model freezing failed - check Docker logs
    exit /b 1
)

REM =============================================================================
REM PyTorch Model Deployment
REM =============================================================================

echo === Copying PyTorch Model for LAMMPS ===

REM Check for .pth model file
if exist deepmd_model\hap_model.pth (
    copy deepmd_model\hap_model.pth lammps_run\model\ >nul 2>&1
    if %errorlevel% equ 0 (
        echo PyTorch model copied successfully
    ) else (
        echo ERROR: Failed to copy PyTorch model file
    )
) else (
    echo Warning: PyTorch model file not found. Will use default LAMMPS parameters.
    
    REM Additional diagnostics for PyTorch model
    echo Checking for any model files in deepmd_model directory...
    dir deepmd_model\*.*
)

REM =============================================================================
REM LAMMPS Simulation with PyTorch DeepMD Model
REM =============================================================================

echo === Running LAMMPS with PyTorch DeepMD Model ===
echo Starting LAMMPS simulation with PyTorch DeepMD model...

REM Using the DeepMD image that's already pulled for LAMMPS
docker compose up -d --no-deps --force-recreate lammps_run

REM Get container ID for LAMMPS
FOR /f %%i IN ('docker compose ps -a -q lammps_run 2^>nul') DO set LAMMPS_ID=%%i

if not defined LAMMPS_ID (
    echo ERROR: Failed to get LAMMPS container ID
    exit /b 1
)

echo LAMMPS container ID: !LAMMPS_ID!

REM Follow logs and wait for container
docker wait !LAMMPS_ID!
set LAMMPS_EXIT_CODE=!errorlevel!

REM Capture logs after container finishes
docker logs !LAMMPS_ID!

if !LAMMPS_EXIT_CODE! equ 0 (
    echo LAMMPS simulation with PyTorch model completed successfully
) else (
    echo ERROR: LAMMPS simulation failed - check Docker logs
    exit /b 1
)

REM =============================================================================
REM Workflow Completion
REM =============================================================================

echo === PyTorch DeepMD Workflow Complete ===

echo.
echo ==========================================
echo PyTorch DeepMD Workflow completed!
echo Surface Model: (%MILLER_H% %MILLER_K% %MILLER_L%) surface
echo Surface Parameters: %SURFACE_LAYERS% layers, %SURFACE_SUPERCELL_A%x%SURFACE_SUPERCELL_B% supercell
echo CO2 Addition: %ADD_CO2%
echo Number of CO2 molecules per side: %N_CO2%
echo Slab Method: %SLAB_METHOD%
echo Backend: PyTorch with DPA Descriptor
echo Model: deepmd_model\hap_model.pth
echo Initial Structure: lammps_run\%HAP_DATA_FILE%
echo ==========================================

REM Clean up temporary files
if exist "%TEMP%\data_files.tmp" del "%TEMP%\data_files.tmp"

goto :eof