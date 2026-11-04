param([ValidatePattern('^[A-Za-z0-9_-]+$')][string]$Run='refine-001')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$runRoot=(Resolve-Path -LiteralPath ".research-runs/$Run").Path
$source=(Resolve-Path -LiteralPath (Join-Path $runRoot 'dft-relax')).Path
$target=[IO.Path]::GetFullPath((Join-Path $runRoot 'dft-relax-attempt3'))
if (-not $source.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar) -or -not $target.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar)) { throw 'Archive path leaves run directory' }
if (Test-Path -LiteralPath $target) { throw 'Preserving prior archive; this comparison runs only once' }
$job=Get-Content -LiteralPath (Join-Path $source 'job.json') -Raw | ConvertFrom-Json
$audit=Get-Content -LiteralPath (Join-Path $runRoot 'dft_relaxation.json') -Raw | ConvertFrom-Json
$marker=Get-Content -LiteralPath (Join-Path $runRoot 'dft-audit.manifest.json') -Raw | ConvertFrom-Json
if ($job.status -ne 'completed' -or $marker.status -ne 'failed' -or -not $audit.force_gate -or -not $audit.chemistry.geometry_gate -or $audit.optimizer_converged) { throw 'Requires reviewed terminal displacement-convergence failure with valid force/chemistry gates' }
$running=& docker ps --filter 'name=hap-refine-' --format '{{.Names}}'
if ($LASTEXITCODE -ne 0 -or $running) { throw 'A refinement calculation remains active or cannot be checked' }
& python workflow/refine_pipeline.py dft-diagnose --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Diagnosis failed; no archive or new calculation launched' }
Move-Item -LiteralPath $source -Destination $target
foreach ($name in @('dft-audit.manifest.json','dft-audit.console.log','dft-audit.console.execution.json')) { Move-Item -LiteralPath (Join-Path $runRoot $name) -Destination (Join-Path $target $name) }
Move-Item -LiteralPath (Join-Path $runRoot 'dft_relaxation.json') -Destination (Join-Path $target 'audit-result.json')
& python workflow/refine_pipeline.py dft-lbfgs-prepare --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'L-BFGS preparation failed; original attempt retained' }
& python workflow/run_refinement.py all --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Bounded optimizer comparison failed; inspect retained results' }
