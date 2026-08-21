#!/usr/bin/env bash
#
# check-sync-fxpower.sh -- verify zzfisher's copy of the fxpower
# shared sources still matches the recorded sync checksums.
#
# This is the copy of tools/check-sync.sh installed INTO zzfisher by
# sync-to-zzfisher.sh (under the "-fxpower" name, since zzfisher
# already carries an unrelated rx2 sync under tools/check-sync.sh
# from rgt47/fisherexacttestrx2). Run from zzfisher's own tools/
# directory.
#
# Exit status 0 = in sync, 1 = drift detected.

set -u

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)

manifest="$script_dir/sync-manifest-fxpower.txt"
checksums="$script_dir/sync-checksums-fxpower.sha256"

if [ ! -f "$manifest" ]; then
  echo "check-sync-fxpower: missing $manifest" >&2
  exit 1
fi

if [ ! -f "$checksums" ]; then
  echo "check-sync-fxpower: missing $checksums" >&2
  echo "check-sync-fxpower: run tools/sync-to-zzfisher.sh in" >&2
  echo "                    rgt47/fisherpowerunequaln to generate it" >&2
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

while read -r line; do
  case "$line" in
    ''|\#*) continue ;;
  esac
  if [ ! -f "$line" ]; then
    echo "check-sync-fxpower: MISSING  $line"
    missing=1
  fi
done < "$manifest"

if [ "$missing" -ne 0 ]; then
  echo "check-sync-fxpower: FAILED (files listed in the manifest are absent)"
  exit 1
fi

while read -r want path; do
  case "$want" in
    ''|\#*) continue ;;
  esac
  got=$(sha_cmd "$path" | awk '{ print $1 }')
  if [ "$got" != "$want" ]; then
    echo "check-sync-fxpower: DRIFT    $path"
    status=1
  fi
done < "$checksums"

n_manifest=$(grep -cvE '^[[:space:]]*(#|$)' "$manifest")
n_checksums=$(grep -cvE '^[[:space:]]*(#|$)' "$checksums")
if [ "$n_manifest" -ne "$n_checksums" ]; then
  echo "check-sync-fxpower: manifest lists $n_manifest files but checksums"
  echo "                    cover $n_checksums; re-run the export"
  status=1
fi

if [ "$status" -eq 0 ]; then
  echo "check-sync-fxpower: OK ($n_manifest shared files match)"
  if [ -f "$script_dir/sync-provenance-fxpower.txt" ]; then
    grep -E '^(source-commit|synced-at):' \
      "$script_dir/sync-provenance-fxpower.txt" || true
  fi
else
  echo "check-sync-fxpower: FAILED -- the two repositories have drifted"
fi

exit "$status"
