#!/usr/bin/env bash
#
# check-sync.sh -- verify that this repository's copy of the shared
# kernel sources still matches the recorded sync checksums.
#
# This script is identical in rgt47/fisherexacttestrx2 and
# rgt47/zzfisher, and each repository verifies only its own files.
# No network access and no checkout of the other repository is
# required: both hold the same tools/sync-checksums.sha256, written
# at the last successful sync.
#
# Exit status 0 = in sync, 1 = drift detected.
#
# A failure here means one of two things:
#
#   1. A shared file was edited in this repository without running
#      the export. If this is the research workspace, run
#      tools/sync-to-zzfisher.sh. If this is the package, the edit
#      belongs upstream: move it to the research workspace instead.
#   2. A shared file was edited in the OTHER repository and the
#      checksums were updated there. Pull, then re-sync.

set -u

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)

manifest="$script_dir/sync-manifest.txt"
checksums="$script_dir/sync-checksums.sha256"

if [ ! -f "$manifest" ]; then
  echo "check-sync: missing $manifest" >&2
  exit 1
fi

if [ ! -f "$checksums" ]; then
  echo "check-sync: missing $checksums" >&2
  echo "check-sync: run tools/sync-to-zzfisher.sh in the research" >&2
  echo "            workspace to generate it" >&2
  exit 1
fi

if command -v sha256sum >/dev/null 2>&1; then
  sha_cmd() { sha256sum "$@"; }
else
  sha_cmd() { shasum -a 256 "$@"; }
fi

cd "$repo_root" || exit 1

status=0
missing=0

# Every manifest path must exist locally.
while read -r line; do
  case "$line" in
    ''|\#*) continue ;;
  esac
  if [ ! -f "$line" ]; then
    echo "check-sync: MISSING  $line"
    missing=1
  fi
done < "$manifest"

if [ "$missing" -ne 0 ]; then
  echo "check-sync: FAILED (files listed in the manifest are absent)"
  exit 1
fi

# Compare recorded checksums against the working tree.
while read -r want path; do
  case "$want" in
    ''|\#*) continue ;;
  esac
  got=$(sha_cmd "$path" | awk '{ print $1 }')
  if [ "$got" != "$want" ]; then
    echo "check-sync: DRIFT    $path"
    status=1
  fi
done < "$checksums"

# Catch files added to the manifest but never checksummed.
n_manifest=$(grep -cvE '^[[:space:]]*(#|$)' "$manifest")
n_checksums=$(grep -cvE '^[[:space:]]*(#|$)' "$checksums")
if [ "$n_manifest" -ne "$n_checksums" ]; then
  echo "check-sync: manifest lists $n_manifest files but checksums"
  echo "            cover $n_checksums; re-run the export"
  status=1
fi

if [ "$status" -eq 0 ]; then
  echo "check-sync: OK ($n_manifest shared files match)"
  if [ -f "$script_dir/sync-provenance.txt" ]; then
    grep -E '^(source-commit|synced-at):' \
      "$script_dir/sync-provenance.txt" || true
  fi
else
  echo "check-sync: FAILED -- the two repositories have drifted"
fi

exit "$status"
