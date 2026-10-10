#!/bin/bash
## Fetch pmsimstats at commit 58b32a9 ('orig'), apply the patches, check that
## each applies without fuzz or rejects, extract the pinned stress-test driver
## used as the reference, and run verify.R. Requires an authenticated GitHub
## CLI. Run from the repository root:
##   bash analysis/patches/orig-to-refined/verify.sh [workdir [--fast]]
## The work directory (default: a new temporary directory) is kept, so the
## patched trees can be inspected afterwards; further arguments are passed to
## verify.R.
set -euo pipefail
here=analysis/patches/orig-to-refined
sha=58b32a92b99a8beb145de00c7bd8a0880f06f0f0
## The stress-test driver as committed before its defaults were changed to
## anchored carryover and randomized allocation: the fixed reference.
ref=32389168e844adb628e1caed8bd62b380828ffec
work=${1:-$(mktemp -d)}
mkdir -p "$work"
echo "work directory: $work"

gh api "repos/rchendrickson/pmsimstats/tarball/$sha" > "$work/orig.tar.gz"
rm -rf "$work/src"
mkdir -p "$work/src"
tar -xzf "$work/orig.tar.gz" -C "$work/src"
src=$(ls -d "$work"/src/*/)

fails=0
mktree() {
  rm -rf "$work/$1"
  mkdir -p "$work/$1"
  cp -R "$src/R" "$src/DESCRIPTION" "$src/data" "$work/$1/"
}
apply_one() {
  local tree=$1 patchfile=$2
  local out gitok
  if (cd "$work/$tree" && git apply --check -p1 "$OLDPWD/$here/$patchfile"); then
    gitok=yes
  else
    gitok=no
  fi
  out=$(cd "$work/$tree" && patch -p1 --no-backup-if-mismatch -i "$OLDPWD/$here/$patchfile" 2>&1) || true
  if [ "$gitok" = yes ] && ! grep -qiE 'fuzz|fail|reject|rej' <<<"$out"; then
    echo "PASS apply $patchfile to $tree (git apply --check clean; patch: no fuzz, no rejects)"
  else
    echo "FAIL apply $patchfile to $tree"
    echo "$out"
    fails=$((fails + 1))
  fi
  grep -i offset <<<"$out" | sed 's/^/     note: /' || true
}

mktree orig
mktree core;     apply_one core orig-to-refined.patch
mktree options;  apply_one options orig-to-refined.patch
                 apply_one options orig-to-refined-options.patch
mktree analysis; apply_one analysis orig-to-refined-analysis.patch
mktree full;     apply_one full orig-to-refined.patch
                 apply_one full orig-to-refined-options.patch
                 apply_one full orig-to-refined-analysis.patch

mkdir -p "$work/ref"
git show "$ref:analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R" \
  > "$work/ref/11-strawman-stress-test.R"
git show "$ref:analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R" \
  > "$work/ref/04-power-simulation.R"
echo "reference driver: commit $ref"

if [ "$fails" -gt 0 ]; then
  echo "FAIL $fails patch application(s); verify.R not run"
  exit 1
fi
Rscript "$here/verify.R" "$work" "${@:2}"
