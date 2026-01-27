#!/bin/bash
set -e

echo "Starting GECKO-A..."
cd "$(dirname "$0")/.."

# Check if Docker daemon is running
if ! docker info >/dev/null 2>&1; then
    echo "❌ Docker daemon is not running."
    echo "Please start Docker Desktop."
    open -a Docker
    # Wait for Docker to start
    echo "Waiting for Docker to start..."
    while ! docker info >/dev/null 2>&1; do
        sleep 2
        echo -n "."
    done
    echo "Docker started."
fi

docker-compose -f docker/docker-compose.yml up -d
echo "Opening browser..."
open "http://localhost:8000"
