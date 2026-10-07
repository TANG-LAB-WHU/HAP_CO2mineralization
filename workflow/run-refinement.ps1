param(
  [ValidateSet('all','references','aimd','dft','fragments','train','lammps','plots')][string]$Stage='all',
  [ValidatePattern('^[A-Za-z0-9_-]+$')][string]$Run='refine-001'
)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
if ($Stage -eq 'plots') {
  & python workflow/plot_refinement.py --run ".research-runs/$Run" --animate
} else {
  & python workflow/run_refinement.py $Stage --run ".research-runs/$Run"
}
if ($LASTEXITCODE -ne 0) { throw "Refinement stage failed ($LASTEXITCODE); evidence/logs retained" }
