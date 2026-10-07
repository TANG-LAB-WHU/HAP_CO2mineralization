param([ValidatePattern('^[A-Za-z0-9_-]+$')][string]$Run='refine-001')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$runRoot=(Resolve-Path -LiteralPath ".research-runs/$Run").Path
$source=(Resolve-Path -LiteralPath (Join-Path $runRoot 'dft-relax')).Path
$target=[IO.Path]::GetFullPath((Join-Path $runRoot 'dft-relax-attempt1'))
if (-not $source.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar) -or -not $target.StartsWith($runRoot+[IO.Path]::DirectorySeparatorChar)) { throw 'Archive path leaves the run directory' }
if (Test-Path -LiteralPath $target) { throw 'Preserving existing archived attempt' }
$job=Get-Content -LiteralPath (Join-Path $source 'job.json') -Raw | ConvertFrom-Json
$execution=Get-Content -LiteralPath (Join-Path $source 'console.execution.json') -Raw | ConvertFrom-Json
if ($job.status -ne 'failed' -or $execution.error -notmatch 'TimeoutExpired') { throw 'This recovery is only for a confirmed terminal wall-time failure' }
$running=& docker ps --filter 'name=hap-refine-' --format '{{.Names}}'
if ($LASTEXITCODE -ne 0 -or $running) { throw 'A refinement container remains active or could not be checked' }
Move-Item -LiteralPath $source -Destination $target
& python workflow/refine_pipeline.py dft-resume --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Restart preparation failed; all prior outputs are retained' }
& python workflow/run_refinement.py all --run ".research-runs/$Run"
if ($LASTEXITCODE -ne 0) { throw 'Continuation failed; inspect retained logs before another action' }
