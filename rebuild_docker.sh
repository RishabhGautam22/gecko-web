#!/bin/bash

# Rebuild Docker Image for GECKO-A
echo "🐳 Building GECKO-A Docker Image..."

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "❌ Docker is not running. Please start Docker Desktop and try again."
    exit 1
fi

# Build the image
# We assume the build context is the project root
docker build -t gecko-a-web -f docker/Dockerfile .

if [ $? -eq 0 ]; then
    echo "✅ Build successful!"
    echo "   You can run it with: docker run -p 8000:8000 -v $(pwd)/data:/data gecko-a-web"
else
    echo "❌ Build failed."
    exit 1
fi
