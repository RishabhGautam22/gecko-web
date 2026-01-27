#!/usr/bin/env bash
# ============================================================================
# GECKO-A Web Interface - Unix Launcher
# Version 3.0.6
# Compatible with: macOS, Linux, WSL
# ============================================================================

set -e

VERSION="3.0.6"
PORT="${PORT:-8000}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Navigate to project root
cd "$(dirname "$0")/.." || exit 1

echo ""
echo -e "${CYAN}========================================================================"
echo "  GECKO-A Web Interface - Unix Launcher"
echo -e "  Version ${VERSION}${NC}"
echo -e "${CYAN}========================================================================${NC}"
echo ""

# Check for Docker
check_docker() {
    if command -v docker &> /dev/null && docker info &> /dev/null; then
        return 0
    fi
    return 1
}

# Docker mode (recommended)
start_docker() {
    echo -e "${GREEN}[OK]${NC} Docker found"
    echo -e "${YELLOW}[INFO]${NC} Starting with Docker (recommended for full functionality)"
    echo ""

    if command -v docker-compose &> /dev/null; then
        docker-compose -f docker/docker-compose.yml up -d
    else
        docker compose -f docker/docker-compose.yml up -d
    fi

    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}[SUCCESS]${NC} GECKO-A is running!"
        echo -e "URL: ${CYAN}http://localhost:${PORT}${NC}"
        echo ""
        echo "To stop: docker-compose -f docker/docker-compose.yml down"

        sleep 2
        # Open browser based on OS
        if [[ "$OSTYPE" == "darwin"* ]]; then
            open "http://localhost:${PORT}"
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            xdg-open "http://localhost:${PORT}" 2>/dev/null || true
        fi
    else
        echo -e "${RED}[ERROR]${NC} Docker failed. Trying Python mode..."
        start_python
    fi
}

# Python mode (limited)
start_python() {
    echo ""
    echo -e "${YELLOW}[WARN]${NC} Starting in LIMITED mode (UI only, no simulations)"
    echo "For full functionality, install Docker"
    echo ""

    # Check Python
    if ! command -v python3 &> /dev/null; then
        echo -e "${RED}[ERROR]${NC} Python 3 not found"
        echo "Please install Python 3.9+ from https://www.python.org"
        exit 1
    fi

    PYVER=$(python3 --version)
    echo -e "${GREEN}[OK]${NC} ${PYVER}"

    # Create venv if needed
    if [ ! -d "venv" ]; then
        echo "Creating virtual environment..."
        python3 -m venv venv
    fi

    # Activate venv
    source venv/bin/activate

    # Install dependencies
    echo "Installing dependencies..."
    pip install -q --upgrade pip
    pip install -q -r gecko_web/requirements.txt

    # Set environment
    export DATA_DIR="$(pwd)/data"
    export GECKO_SOURCE_DIR="$(pwd)/docker/gecko_source"
    export BOXMODEL_SOURCE_DIR="$(pwd)/docker/boxmodel_source"

    # Create data directory
    mkdir -p data/output

    # Kill existing process on port
    if lsof -ti :${PORT} &>/dev/null; then
        echo "Freeing port ${PORT}..."
        lsof -ti :${PORT} | xargs kill -9 2>/dev/null || true
        sleep 1
    fi

    echo ""
    echo -e "${GREEN}[SUCCESS]${NC} Starting GECKO-A Web Interface"
    echo -e "URL: ${CYAN}http://localhost:${PORT}${NC}"
    echo "Press Ctrl+C to stop"
    echo ""

    # Open browser
    sleep 2 && {
        if [[ "$OSTYPE" == "darwin"* ]]; then
            open "http://localhost:${PORT}"
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            xdg-open "http://localhost:${PORT}" 2>/dev/null || true
        fi
    } &

    uvicorn gecko_web.main:app --host 127.0.0.1 --port ${PORT}
}

# Main
if check_docker; then
    start_docker
else
    echo -e "${YELLOW}[WARN]${NC} Docker not available"
    start_python
fi
