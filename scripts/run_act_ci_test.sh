#!/usr/bin/env bash
set -euo pipefail

WORKFLOW="${WORKFLOW:-.github/workflows/ci.yml}"
JOB="${JOB:-test}"
IMAGE="${IMAGE:-catthehacker/ubuntu:act-latest}"
ARCH="${ARCH:-linux/amd64}"
TMP_ROOT="${TMP_ROOT:-${HOME}/.mpi-tmp}"

echo "Cleaning old act containers..."
docker rm -f $(docker ps -aq --filter "name=act-") 2>/dev/null || true

echo "Cleaning Docker cache..."
docker system prune -af
docker volume prune -f

echo "Cleaning act cache..."
rm -rf "${HOME}/.cache/act"
rm -rf "${TMP_ROOT}"/act* 2>/dev/null || true

echo "Disk usage after cleanup:"
df -h
docker system df

echo "Running act workflow..."
echo "  workflow: ${WORKFLOW}"
echo "  job:      ${JOB}"
echo "  image:    ${IMAGE}"
echo "  arch:     ${ARCH}"

act push \
  -W "${WORKFLOW}" \
  -j "${JOB}" \
  --container-architecture "${ARCH}" \
  --bind \
  -P "ubuntu-latest=${IMAGE}"