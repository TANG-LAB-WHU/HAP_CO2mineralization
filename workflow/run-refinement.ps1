param(
  [ValidateSet('all','references','aimd','dft','fragments','train','lammps','plots')][string]$Stage='all',
  [ValidatePattern('^[A-Za-z0-9_-]+$')][string]$Run='refine-001'
)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
if ($Stage -eq 'plots') {
  & python workflow/plot_refinement.py --run ".research-runs/$Run" --animate
} else {
  $active = @(& docker ps --filter 'name=hap-refine-' --format '{{.Names}}')
  if ($LASTEXITCODE -ne 0) { throw 'Cannot inspect Docker jobs; refusing to start a calculation' }
  if ($active.Count -gt 0) { throw "Refinement calculation already active: $($active -join ', '). Inspect that runner instead of starting a duplicate." }
  & python workflow/run_refinement.py $Stage --run ".research-runs/$Run"
}
if ($LASTEXITCODE -ne 0) { throw "Refinement stage failed ($LASTEXITCODE); evidence/logs retained" }
