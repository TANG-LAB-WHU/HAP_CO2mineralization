#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INSTALL=false
WITH_PDF=false

for arg in "$@"; do
  case "$arg" in
    --install) INSTALL=true ;;
    --with-pdf) WITH_PDF=true ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done

cd "$REPO_ROOT"

if [[ ! -d .git ]]; then
  echo "Run this script from inside the cloned HAP_CO2mineralization repository." >&2
  exit 1
fi

if [[ -f research-workflow.gitignore ]]; then
  touch .gitignore
  while IFS= read -r rule || [[ -n "$rule" ]]; do
    [[ -z "$rule" ]] && continue
    grep -Fqx "$rule" .gitignore || printf '%s\n' "$rule" >> .gitignore
  done < research-workflow.gitignore
fi

if ! command -v brew >/dev/null 2>&1; then
  echo "Homebrew is required for the automated install path: https://brew.sh/" >&2
  echo "Install Homebrew, then rerun this script with --install." >&2
  exit 1
fi

if [[ "$INSTALL" == true ]]; then
  command -v python3 >/dev/null 2>&1 || brew install python
  command -v git-lfs >/dev/null 2>&1 || brew install git-lfs
  command -v quarto >/dev/null 2>&1 || brew install --cask quarto
fi

missing=()
for command_name in git python3 quarto; do
  command -v "$command_name" >/dev/null 2>&1 || missing+=("$command_name")
done

if (( ${#missing[@]} > 0 )); then
  echo "Missing commands: ${missing[*]}" >&2
  echo "Rerun with --install after Homebrew is available." >&2
  exit 1
fi

if command -v git-lfs >/dev/null 2>&1; then
  git lfs install --local
fi

if [[ "$WITH_PDF" == true ]]; then
  quarto install tinytex --no-prompt
fi

python3 scripts/check_manuscript.py
quarto check

echo
echo "Local research workflow is ready."
echo "Next: make paper"
echo "Live preview: make preview"

