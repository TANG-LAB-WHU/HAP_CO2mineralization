param([ValidatePattern('^[A-Za-z0-9_-]+$')][string]$Run='refine-001')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$runRoot=(Resolve-Path -LiteralPath ".research-runs/$Run").Path
$source=(Resolve-Path -LiteralPath (Join-Path $runRoot 'dft-relax')).Path
$target=[IO.Path]::GetFullPath((Join-Path $runRoot 'dft-relax-attempt2'))
if (-not $source.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar) -or -not $target.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar)) { throw 'Archive path leaves the run directory' }
if (Test-Path -LiteralPath $target) { throw 'Preserving existing archive; extension can run only once' }
$job=Get-Content -LiteralPath (Join-Path $source 'job.json') -Raw | ConvertFrom-Json
$log=Get-Content -LiteralPath (Join-Path $source 'opt.log') -Raw
$audit=Get-Content -LiteralPath (Join-Path $runRoot 'dft-audit.manifest.json') -Raw | ConvertFrom-Json
if ($job.status -ne 'completed' -or $audit.status -ne 'failed' -or $log -notmatch 'MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED' -or $log -match 'GEOMETRY OPTIMIZATION COMPLETED') { throw 'Requires a terminal step-limit failure and retained failed audit' }
$running=& docker ps --filter 'name=hap-refine-' --format '{{.Names}}'
if ($LASTEXITCODE -ne 0 -or $running) { throw 'Refinement process remains active or cannot be checked' }
Move-Item -LiteralPath $source -Destination $target
foreach ($name in @('dft-audit.manifest.json','dft-audit.console.log','dft-audit.console.execution.json')) { Move-Item -LiteralPath (Join-Path $runRoot $name) -Destination (Join-Path $target $name) }
Move-Item -LiteralPath (Join-Path $runRoot 'dft_relaxation.json') -Destination (Join-Path $target 'audit-result.json')
& python workflow/refine_pipeline.py dft-extend --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Bounded extension preparation failed; archived data retained' }
& python workflow/run_refinement.py all --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Bounded continuation failed; inspect retained records' }
