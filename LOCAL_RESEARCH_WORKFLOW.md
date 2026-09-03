# Local research workflow on macOS

This layer adds a reproducible writing and evidence workflow without changing the existing computation pipeline. HPC submission remains disabled.

## First setup

```bash
brew install git-lfs
git lfs install
mkdir -p ~/Projects
cd ~/Projects
git clone --branch shawn_dev --single-branch \
  https://github.com/TANG-LAB-WHU/HAP_CO2mineralization.git
cd HAP_CO2mineralization
bash scripts/bootstrap_macos.sh --install --with-pdf
code .
```

If the repository is already cloned:

```bash
cd ~/Projects/HAP_CO2mineralization
git fetch origin
git switch shawn_dev
git pull --ff-only origin shawn_dev
bash scripts/bootstrap_macos.sh --install --with-pdf
```

## Daily loop

1. Update `research/research_state.yml` with the current question and next decision.
2. Log searched and screened sources in `literature/search_log.md`.
3. Add verified metadata to `paper/references.bib`.
4. Connect each manuscript claim to a source or result in `literature/evidence.csv`.
5. Edit `paper/manuscript.qmd`.
6. Run `make check` and `make paper`.
7. Review `git diff`, commit, and push `shawn_dev`.

```bash
make check
make paper
git status
git diff
git add AGENTS.md LOCAL_RESEARCH_WORKFLOW.md Makefile research-workflow.gitignore \
  paper literature research scripts workflow .vscode .github .gitignore
git commit -m "feat: add local research writing workflow"
git push origin shawn_dev
```

Open the generated HTML at `paper/_output/manuscript.html`. Use `make preview` for live preview and `make paper-pdf` for PDF.

## Codex role

Install the official Codex extension in VS Code, sign in, open this repository, and keep changes local and reviewable. `AGENTS.md` gives Codex the repository-specific rules.

A good first prompt is:

> Read AGENTS.md, research/research_state.yml, literature/evidence.csv, paper/manuscript.qmd, and the existing project README. Propose one bounded improvement to the manuscript. Do not invent results or citations. Mark missing evidence explicitly, run make check, and show me the diff before committing.

Codex is the implementation and synthesis layer, not the source of evidence. Literature metadata and scientific claims still need verification.

## Future compute handoff

`workflow/profiles/whu-hpc.example.yaml` reserves the Slurm interface but is disabled and contains no credentials. Later, a machine-local untracked file named `workflow/profiles/whu-hpc.local.yaml` can hold non-secret host aliases and allocation settings. Credentials remain in EasyConnect, macOS Keychain, and `~/.ssh/`, never in Git.
