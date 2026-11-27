#!/usr/bin/env bash
# Run before committing/pushing: fails if tracked or staged files look like personal genomic data.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
bad=0
files=$( { git ls-files; git diff --cached --name-only; } | sort -u)
while IFS= read -r f; do
  [ -f "$f" ] || continue
  case "$f" in
    *.vcf|*.vcf.gz|*.cram|*.crai|*.bam|*.bai|*.pgen|*.pvar|*.psam|*.bed|*.bim|*.fam|*.sscore|*.json|*.xlsx|*.pkl|*.fa|*.fasta)
      echo "forbidden file type: $f"; bad=1;;
  esac
  if [ "$(wc -c < "$f" | tr -d " ")" -gt 2000000 ]; then echo "file larger than 2 MB: $f"; bad=1; fi
done <<< "$files"
# Strings that should never appear in public code (edit the patterns for your own identifiers, locally)
if [ -f .personal_patterns ]; then
  if git grep -n -I -f .personal_patterns -- . ':!.personal_patterns' >/dev/null 2>&1; then
    git grep -n -I -f .personal_patterns -- . ':!.personal_patterns' | head; bad=1
  fi
fi
[ $bad -eq 0 ] && echo "OK: no personal-data patterns found" || { echo "FAILED"; exit 1; }
