#!/bin/bash
set -e

echo "Building test environment..."
cd "$(dirname "$0")"
docker-compose -f docker/docker-compose.yml build

echo "Running tests inside container..."
docker-compose -f docker/docker-compose.yml run --rm gecko-app pytest tests/

echo "Tests completed."
