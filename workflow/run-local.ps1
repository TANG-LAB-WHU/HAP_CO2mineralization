param(
  [ValidateSet('prepare','mace','cp2k','check','environment')][string]$Stage='check',
  [string]$Run='smoke-001', [int]$Steps=12
)
$ErrorActionPreference='Stop'
$repo=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repo
if ($Run -notmatch '^[A-Za-z0-9_-]+$') { throw 'Run must be a simple directory name' }
$deepmd='deepmodeling/deepmd-kit@sha256:23d23b922459eb49b23f9d4977b6e453fd17feb2ee03a280c70fb7f4a6452e8c'
$cp2k='cp2k/cp2k@sha256:e2912cc9500a42a08b3e057e383182b96f395b7477c8eb55e93db6e29cb1c2ae'
$common=@('run','--rm','--pull','never','--cpus','4','--memory','8g','--mount',"type=bind,source=$repo,target=/work",'-w','/work','-e','OMP_NUM_THREADS=4')
function Invoke-CheckedDocker([string[]]$Arguments) {
  & docker @Arguments
  if ($LASTEXITCODE -ne 0) { throw "Docker failed: $LASTEXITCODE" }
}
if ($Stage -eq 'prepare') {
  if (Test-Path ".research-runs/$Run/prepare.manifest.json") { throw 'Choose a fresh run name; preserving existing evidence' }
  Invoke-CheckedDocker ($common+@('--entrypoint','/work/.research-runs/venv/bin/python',$deepmd,'workflow/local_compute.py','prepare','--out',".research-runs/$Run"))
} elseif ($Stage -eq 'mace') {
  if (Test-Path ".research-runs/$Run/mace.manifest.json") { throw 'Choose a fresh run name; preserving existing evidence' }
  Invoke-CheckedDocker ($common+@('--gpus','device=0','--entrypoint','/work/.research-runs/venv/bin/python',$deepmd,'workflow/local_compute.py','mace','--out',".research-runs/$Run",'--model','.research-runs/models/mace-mp-small.model','--steps',"$Steps"))
} elseif ($Stage -eq 'cp2k') {
  Invoke-CheckedDocker ($common+@('--entrypoint','/work/.research-runs/venv/bin/python',$deepmd,'workflow/local_compute.py','cp2k-inputs','--out',".research-runs/$Run"))
  foreach ($system in @('co2','slab','adsorbed')) {
    if (Test-Path ".research-runs/$Run/$system.cp2k.log") { throw "Preserving $system log; choose fresh run or explicitly inspect/restart it" }
    Invoke-CheckedDocker ($common+@('--gpus','device=0','-w',"/work/.research-runs/$Run",$cp2k,'cp2k','-i',"$system.inp",'-o',"$system.cp2k.log"))
  }
  & python workflow/cp2k_validate.py ".research-runs/$Run"
  if ($LASTEXITCODE -ne 0) { throw 'CP2K validation failed' }
 } elseif ($Stage -eq 'check') {
  # Existing evidence checks intentionally use POSIX directory-descriptor safeguards.
  Invoke-CheckedDocker ($common+@('--entrypoint','/bin/bash',$deepmd,'-lc','set -e; /work/.research-runs/venv/bin/python scripts/check_manuscript.py; /work/.research-runs/venv/bin/python -m unittest discover -s tests -p "test_*.py"'))
} else {
  Invoke-CheckedDocker ($common+@('--entrypoint','/bin/bash',$deepmd,'-lc','/work/.research-runs/venv/bin/python -m pip freeze > .research-runs/environment/mace-freeze.txt'))
  & docker image inspect $deepmd $cp2k | Set-Content .research-runs/environment/images.json
  & nvidia-smi -q | Set-Content .research-runs/environment/nvidia-smi.txt
  & git rev-parse HEAD | Set-Content .research-runs/environment/source-commit.txt
  & git diff | Set-Content .research-runs/environment/source.patch
}
