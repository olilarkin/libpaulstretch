#!/usr/bin/env bash
#
# Cut a release: bump npm/package.json, commit, tag, push.
#
# Usage:
#   scripts/release.sh patch          # 0.2.0 -> 0.2.1
#   scripts/release.sh minor          # 0.2.0 -> 0.3.0
#   scripts/release.sh major          # 0.2.0 -> 1.0.0
#   scripts/release.sh 0.4.2          # explicit version
#
# After push, two workflows fire on the v* tag:
#   - Release       (npm/wasm to GitHub Packages)
#   - Release CLI   (paulstretch binaries to a GitHub Release)
#
# The script aborts on a dirty tree, a non-main branch, an out-of-date
# remote, or a duplicate tag. It pauses for confirmation before the final
# push so you can sanity-check the commit and tag locally first.

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <patch|minor|major|X.Y.Z>" >&2
    exit 1
fi

ARG="$1"
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

# --- Preflight ----------------------------------------------------------------

if [ -n "$(git status --porcelain)" ]; then
    echo "error: working tree is not clean:" >&2
    git status --short >&2
    exit 1
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [ "$BRANCH" != "main" ]; then
    echo "error: must be on main (currently on '$BRANCH')" >&2
    exit 1
fi

echo "Fetching origin..."
git fetch --quiet origin

LOCAL="$(git rev-parse @)"
REMOTE="$(git rev-parse '@{u}')"
if [ "$LOCAL" != "$REMOTE" ]; then
    echo "error: local main is not in sync with origin/main" >&2
    echo "  local:  $LOCAL" >&2
    echo "  origin: $REMOTE" >&2
    exit 1
fi

# --- Compute new version ------------------------------------------------------

case "$ARG" in
    patch|minor|major)
        # npm version prints "vX.Y.Z" on its own line; strip the leading v.
        NEW_VERSION="$(cd npm && npm version "$ARG" --no-git-tag-version | tail -1)"
        NEW_VERSION="${NEW_VERSION#v}"
        ;;
    *)
        if ! [[ "$ARG" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
            echo "error: '$ARG' is not patch|minor|major or a semver X.Y.Z" >&2
            exit 1
        fi
        (cd npm && npm version "$ARG" --no-git-tag-version --allow-same-version >/dev/null)
        NEW_VERSION="$ARG"
        ;;
esac

TAG="v$NEW_VERSION"

# Reject duplicate tags before we make any local commits.
if git rev-parse --verify --quiet "refs/tags/$TAG" >/dev/null; then
    echo "error: tag $TAG already exists locally" >&2
    git checkout -- npm/package.json
    exit 1
fi
if git ls-remote --tags origin "refs/tags/$TAG" | grep -q "$TAG"; then
    echo "error: tag $TAG already exists on origin" >&2
    git checkout -- npm/package.json
    exit 1
fi

# --- Commit + tag -------------------------------------------------------------

git add npm/package.json
git commit -m "Bump version to $TAG"
git tag -a "$TAG" -m "$TAG"

echo
echo "Ready to push:"
echo "  commit: $(git rev-parse --short HEAD) — $(git log -1 --format=%s)"
echo "  tag:    $TAG"
echo

read -r -p "Push to origin? [y/N] " ANSWER
case "$ANSWER" in
    y|Y|yes)
        git push origin main
        git push origin "$TAG"
        echo
        echo "Pushed. Watch the release workflows:"
        echo "  gh run list --limit 5"
        ;;
    *)
        echo
        echo "Aborted before push. To undo locally:"
        echo "  git tag -d $TAG"
        echo "  git reset --hard HEAD~1"
        exit 1
        ;;
esac
