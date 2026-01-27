#!/bin/bash
set -e

echo "=== GECKO-A Installer for macOS ==="

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# 1. Check for Docker
if command_exists docker; then
    echo "✅ Docker is installed."
else
    echo "❌ Docker is not installed."
    echo "Attempting to install Docker via Homebrew..."
    
    if ! command_exists brew; then
        echo "Homebrew is not installed. Please install Docker Desktop manually from: https://www.docker.com/products/docker-desktop/"
        exit 1
    fi

    brew install --cask docker
    echo "Docker installed. Please open Docker Desktop from your Applications folder to start the Docker engine."
    read -p "Press Enter once Docker Desktop is running..."
fi

# Check if Docker daemon is running
if ! docker info >/dev/null 2>&1; then
    echo "❌ Docker daemon is not running."
    echo "Please start Docker Desktop and try again."
    exit 1
fi

# 2. Build the application
echo "Building GECKO-A application..."
cd "$(dirname "$0")/.."
docker-compose -f docker/docker-compose.yml build

# 3. Start the application
echo "Starting GECKO-A application..."
docker-compose -f docker/docker-compose.yml up -d

echo "✅ Installation complete!"
echo "The application is running at http://localhost:8000"
open "http://localhost:8000"
