@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM Four-stage unified workflow for 2.ML_AIMD
REM 1) AIMD (CP2K)
REM 2) Training dataset synthesis
REM 3) DeepMD training + freeze
REM 4) LAMMPS large-scale simulation
REM ============================================================

set "HKL_PARAM_TRAIN=112"
set "HKL_PARAM_LAMMPS=%HKL_PARAM_TRAIN%"
if not "%~1"=="" set "HKL_PARAM_TRAIN=%~1"
if not "%~2"=="" set "HKL_PARAM_LAMMPS=%~2"

set "HAP_PREFIX=HAP_%HKL_PARAM_TRAIN%"
set "STEP1_DIR=Step1_aimd_cp2k_runs\%HAP_PREFIX%_Perfect"
set "STEP2_DIR=Step2_dataset_synthesis"
set "STEP3_DIR=Step3_mlip_deepmd"
set "STEP4_DIR=Step4_lammps_scaleup"

set "POOL_DIR=%STEP2_DIR%\deepmd_pool"
set "POOL_TARGET=%POOL_DIR%\%HAP_PREFIX%"
set "MULTI_INPUT=%STEP3_DIR%\input.multisystem.json"
set "PYTHON_CMD="

echo ==========================================
echo 2.ML_AIMD Unified Workflow
echo Training HKL: %HKL_PARAM_TRAIN%
echo LAMMPS HKL: %HKL_PARAM_LAMMPS%
echo Training structure prefix: %HAP_PREFIX%
echo ==========================================

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
    echo ERROR: Python not found in PATH.
    exit /b 1
)

if not exist "%STEP1_DIR%" (
    echo ERROR: Step1 directory not found: %STEP1_DIR%
    exit /b 1
)
if not exist "%STEP2_DIR%\generate_deepmd_data.py" (
    echo ERROR: Missing Step2 script: %STEP2_DIR%\generate_deepmd_data.py
    exit /b 1
)
if not exist "%STEP3_DIR%\prepare_multisystem_input.py" (
    echo ERROR: Missing Step3 script: %STEP3_DIR%\prepare_multisystem_input.py
    exit /b 1
)
if not exist "%STEP3_DIR%\input_template.json" (
    echo ERROR: Missing Step3 template: %STEP3_DIR%\input_template.json
    exit /b 1
)
if not exist "%STEP4_DIR%\in.deepmd.lammps" (
    echo ERROR: Missing Step4 input: %STEP4_DIR%\in.deepmd.lammps
    exit /b 1
)

REM ============================================================
REM Stage 1: AIMD (CP2K)
REM ============================================================
echo === Stage 1/4: AIMD (CP2K) ===

set "OPT_TRAJ=%STEP1_DIR%\%HAP_PREFIX%_opt-pos-1.xyz"
set "FINAL_COORD=%STEP1_DIR%\%HAP_PREFIX%_final.xyz"
set "EXTRACT_SCRIPT=%STEP1_DIR%\extract_last_frame.py"

REM Build final coordinate file from optimization trajectory if needed
if exist "%OPT_TRAJ%" (
    if not exist "%EXTRACT_SCRIPT%" (
        echo ERROR: Missing extractor script: %EXTRACT_SCRIPT%
        exit /b 1
    )
    if not exist "%FINAL_COORD%" (
        echo Building final coordinate from optimization trajectory...
        "%PYTHON_CMD%" "%EXTRACT_SCRIPT%" "%OPT_TRAJ%" "%FINAL_COORD%"
        if %errorlevel% neq 0 (
            echo ERROR: Failed to build %FINAL_COORD%
            exit /b 1
        )
    )
) else (
    if not exist "%FINAL_COORD%" (
        echo ERROR: Missing both optimization trajectory and final coordinates:
        echo   %OPT_TRAJ%
        echo   %FINAL_COORD%
        exit /b 1
    )
)

if not exist "%STEP1_DIR%\%HAP_PREFIX%_md-pos.xyz" (
    echo Running CP2K for %HAP_PREFIX% ...
    set "HKL_PARAM=%HKL_PARAM_TRAIN%"
    docker compose up -d --no-deps --force-recreate cp2k_run
    for /f %%i in ('docker compose ps -a -q cp2k_run 2^>nul') do set CP2K_ID=%%i
    if not defined CP2K_ID (
        echo ERROR: Could not get cp2k_run container id
        exit /b 1
    )
    docker wait !CP2K_ID!
    set CP2K_EXIT=!errorlevel!
    docker logs !CP2K_ID!
    if !CP2K_EXIT! neq 0 (
        echo ERROR: CP2K stage failed.
        exit /b 1
    )
) else (
    echo CP2K trajectory already exists. Skip CP2K run.
)

if not exist "%STEP1_DIR%\%HAP_PREFIX%_md-pos.xyz" (
    echo ERROR: Missing CP2K trajectory after Stage 1: %STEP1_DIR%\%HAP_PREFIX%_md-pos.xyz
    exit /b 1
)

REM ============================================================
REM Stage 2: Synthesize training dataset
REM ============================================================
echo === Stage 2/4: Synthesizing DeepMD datasets ===

"%PYTHON_CMD%" "%STEP2_DIR%\generate_deepmd_data.py" --cp2k-source-dir "%STEP1_DIR%" --output-dir "%POOL_TARGET%" --project-name "%HAP_PREFIX%"
if %errorlevel% neq 0 (
    echo ERROR: Failed to generate DeepMD dataset for %HAP_PREFIX%.
    exit /b 1
)

"%PYTHON_CMD%" "%STEP3_DIR%\prepare_multisystem_input.py" --pool-dir "%POOL_DIR%" --output "%MULTI_INPUT%" --template "%STEP3_DIR%\input_template.json" --container-prefix "/mnt/deepmd_pool"
if %errorlevel% neq 0 (
    echo ERROR: Failed to prepare multi-system training input.
    exit /b 1
)

REM ============================================================
REM Stage 3: DeepMD training and freeze
REM ============================================================
echo === Stage 3/4: Training DeepMD and exporting model ===

set "DEEPMD_INPUT_FILE=input.multisystem.json"
docker compose up -d --no-deps --force-recreate deepmd_train
for /f %%i in ('docker compose ps -a -q deepmd_train 2^>nul') do set TRAIN_ID=%%i
if not defined TRAIN_ID (
    echo ERROR: Could not get deepmd_train container id
    exit /b 1
)
docker wait !TRAIN_ID!
set TRAIN_EXIT=!errorlevel!
docker logs !TRAIN_ID!
if !TRAIN_EXIT! neq 0 (
    echo ERROR: DeepMD training failed.
    exit /b 1
)

docker compose up -d --no-deps --force-recreate deepmd_freeze
for /f %%i in ('docker compose ps -a -q deepmd_freeze 2^>nul') do set FREEZE_ID=%%i
if not defined FREEZE_ID (
    echo ERROR: Could not get deepmd_freeze container id
    exit /b 1
)
docker wait !FREEZE_ID!
set FREEZE_EXIT=!errorlevel!
docker logs !FREEZE_ID!
if !FREEZE_EXIT! neq 0 (
    echo ERROR: DeepMD freeze failed.
    exit /b 1
)

if not exist "%STEP3_DIR%\model\hap_model.pth" (
    echo ERROR: Expected model not found: %STEP3_DIR%\model\hap_model.pth
    exit /b 1
)

REM ============================================================
REM Stage 4: LAMMPS large-scale simulation
REM ============================================================
echo === Stage 4/4: LAMMPS scale-up simulation ===

if not exist "%STEP4_DIR%\create_hap_slabs.py" (
    echo ERROR: Missing slab generator: %STEP4_DIR%\create_hap_slabs.py
    exit /b 1
)

pushd "%STEP4_DIR%"
"%PYTHON_CMD%" create_hap_slabs.py %HKL_PARAM_LAMMPS:~0,1% %HKL_PARAM_LAMMPS:~1,1% %HKL_PARAM_LAMMPS:~2,1% --layers 4 --vacuum 15 --supercell 3 3 --method ase --add-co2 --n-co2 1 --co2-height 10 --output hap_%HKL_PARAM_LAMMPS%
if %errorlevel% neq 0 (
    popd
    echo ERROR: Failed to generate LAMMPS slab data.
    exit /b 1
)
popd

docker compose up -d --no-deps --force-recreate lammps_run
for /f %%i in ('docker compose ps -a -q lammps_run 2^>nul') do set LAMMPS_ID=%%i
if not defined LAMMPS_ID (
    echo ERROR: Could not get lammps_run container id
    exit /b 1
)
docker wait !LAMMPS_ID!
set LAMMPS_EXIT=!errorlevel!
docker logs !LAMMPS_ID!
if !LAMMPS_EXIT! neq 0 (
    echo ERROR: LAMMPS stage failed.
    exit /b 1
)

echo ==========================================
echo Workflow completed successfully for %HAP_PREFIX%
echo LAMMPS HKL used for scale-up: %HKL_PARAM_LAMMPS%
echo Stage1: %STEP1_DIR%
echo Stage2 pooled system: %POOL_TARGET%
echo Stage3 model: %STEP3_DIR%\model\hap_model.pth
echo Stage4 outputs: %STEP4_DIR%
echo ==========================================