# GECKO-A Web Interface - Deployment & Docker Instructions
## Atmospheric Chemistry Mechanism Generator & Box Model Simulator

**Application Version:** 3.0.6
**Primary Stack:** FastAPI (Python 3.9+) + Fortran 90 (GECKO-A & BOXMODEL4GECKO)
**Author:** Deeksha Sharma

---

## Quick Start - Fly.io Deployment (Recommended)

**For detailed Fly.io instructions, see: [FLY_DEPLOYMENT.md](FLY_DEPLOYMENT.md)**

```bash
# 1. Install Fly CLI
curl -L https://fly.io/install.sh | sh

# 2. Login
fly auth login

# 3. Launch app
fly launch --no-deploy

# 4. Create storage volume
fly volumes create gecko_data --region iad --size 10

# 5. Deploy (takes 45-90 min first time)
fly deploy

# 6. Open
fly open
```

**Live URL:** `https://your-app-name.fly.dev`

---

## Table of Contents
1. [Application Analysis Report](#phase-1-application-analysis-report)
2. [Docker Optimization](#phase-2-docker-optimization)
3. [Cross-Platform Compatibility](#phase-3-cross-platform-compatibility)
4. [Dependency Management](#phase-4-dependency-management)
5. [GitHub Repository Setup](#phase-5-github-repository-setup)
6. [CI/CD Workflows](#phase-6-cicd-workflows)
7. [Cloud Deployment](#phase-7-cloud-deployment)
8. [Execution Checklist](#phase-8-execution-checklist)

---

## Phase 1: Application Analysis Report

### 1.1 Project Overview

```
APPLICATION ANALYSIS REPORT
===========================
Project Name: GECKO-A Web Interface
Primary Language: Python 3.9+ (FastAPI) + Fortran 90
Framework(s): FastAPI, Jinja2, RDKit, 3Dmol.js, Cytoscape.js
Package Manager: pip (requirements.txt)
Entry Point: gecko_web/main.py
Default Port: 8000
Database: None (in-memory compound database, JSON job persistence)
External Services: CDN (3Dmol.js, Cytoscape.js, SMILES Drawer)

FORTRAN COMPONENTS:
- GECKO-A: Atmospheric mechanism generator
- BOXMODEL4GECKO: Kinetic box model solver
- Both require: gfortran, NetCDF-Fortran libraries

DEPENDENCIES SUMMARY:
- Python packages: 38 total
- Scientific: NumPy, SciPy, Pandas, RDKit, NetworkX
- Web: FastAPI, Uvicorn, Jinja2
- Visualization: Matplotlib, Graphviz, Pillow, ReportLab
- Data: NetCDF4, Xarray
- Platform-specific: RDKit (binary wheels), NetCDF libraries (system)

POTENTIAL ISSUES:
- Fortran filename buffer overflow (100 char limit) - mitigated via /tmp/gk/
- RDKit requires specific platform wheels
- Graphviz requires system binary installation
- NetCDF requires system library installation
- Large build context (~1.3GB with venv/Windows installer)
```

### 1.2 Directory Structure

```
GECKO/
├── gecko_web/                    # Main Python application
│   ├── main.py                   # FastAPI entry point (3,044 lines)
│   ├── mechanism_diagram.py      # Mechanism parsing & diagrams
│   ├── pathway_visualizer.py     # RDKit + Graphviz visualization
│   ├── enhanced_visualizer.py    # Advanced color visualization
│   ├── reaction_tree.py          # Reaction pathway analysis
│   ├── postprocessing.py         # Dynamic partitioning
│   ├── mass_balance.py           # Atom conservation
│   ├── combined_workflow.py      # Generator + Box Model
│   ├── requirements.txt          # Python dependencies
│   ├── chemdata/                 # Chemical database
│   │   ├── compound_database.py  # 150+ VOC compounds
│   │   ├── voc_categories.py     # VOC categorization
│   │   └── reaction_data.py      # Kinetics data
│   ├── templates/
│   │   └── index.html            # SPA template
│   └── static/
│       ├── css/styles.css
│       └── js/
│           ├── app.js            # Main UI (48.6 KB)
│           ├── api.js            # REST client
│           ├── reaction-tree.js  # Cytoscape integration
│           └── structure3d.js    # 3Dmol integration
│
├── docker/                       # Docker & Fortran sources
│   ├── Dockerfile
│   ├── docker-compose.yml
│   ├── docker_arch.gnu           # NetCDF architecture config
│   ├── gecko_source/             # GECKO-A Fortran (80+ files)
│   │   ├── LIB/                  # Fortran 90 modules
│   │   ├── DATA/                 # Reaction data
│   │   ├── INPUT/                # gecko.nml, cheminput.dat
│   │   └── OBJ/makefile
│   └── boxmodel_source/          # BOXMODEL4GECKO Fortran
│       ├── LIBSRC/               # Fortran subroutines
│       ├── PROG/                 # Main executables
│       ├── CHEMDAT/              # Chemical mechanisms (226 dirs)
│       └── build.sh
│
├── tests/                        # Pytest test suite
├── data/                         # Job outputs (mounted volume)
├── samples/                      # Sample input files
├── scripts/                      # Utility scripts
│
├── start_gecko.command           # macOS launcher
├── rebuild_docker.sh             # Docker build script
├── compile_local.sh              # Local Fortran compilation
├── run_tests.sh                  # Test runner
│
├── README.md
├── DOCUMENTATION.md
└── CHANGELOG.md
```

### 1.3 Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `DATA_DIR` | Job output directory | `/data` (Docker) or `./data` (local) |
| `GECKO_SOURCE_DIR` | GECKO-A source location | `/app/gecko_source` |
| `BOXMODEL_SOURCE_DIR` | Box Model source | `/app/boxmodel_source` |
| `PORT` | Application port | `8000` |

### 1.4 API Endpoints Summary

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/` | GET | Web UI (index.html) |
| `/api/jobs` | GET/POST | Job management |
| `/api/jobs/{id}` | GET/DELETE | Job status/removal |
| `/api/jobs/{id}/results` | GET | Job results |
| `/api/compounds` | GET | Compound database |
| `/api/compounds/search/{q}` | GET | Search compounds |
| `/api/kinetics/{name}` | GET | Rate constants |
| `/api/workflow/combined` | POST | Generator + Box Model |
| `/api/environment` | GET | System status |
| `/health` | GET | Health check (needs creation) |

---

## Phase 2: Docker Optimization

### 2.1 Current Dockerfile Issues

1. **No multi-stage build** - Single stage increases image size
2. **No .dockerignore** - Copies unnecessary files (venv, Windows installer = 1GB+)
3. **No non-root user** - Security concern
4. **No health check** - Container orchestration limitation
5. **No layer caching optimization** - Slow rebuilds
6. **Missing Graphviz system binary** - Required for diagrams

### 2.2 Optimized Dockerfile

Create `docker/Dockerfile.optimized`:

```dockerfile
# =============================================================================
# GECKO-A Web Interface - Optimized Multi-Stage Dockerfile
# Version: 3.0.6
# Platforms: linux/amd64, linux/arm64
# =============================================================================

# -----------------------------------------------------------------------------
# Stage 1: Fortran Build Stage
# -----------------------------------------------------------------------------
FROM python:3.9-slim AS fortran-builder

# Install Fortran build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    make \
    gawk \
    libnetcdf-dev \
    libnetcdff-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy and compile GECKO-A source
COPY docker/gecko_source ./gecko_source
RUN cd gecko_source/OBJ && make clean && make \
    && mkdir -p ../RUN/OUT

# Copy and compile BOXMODEL source
COPY docker/boxmodel_source ./boxmodel_source
COPY docker/docker_arch.gnu ./boxmodel_source/SCRIPTS/arch/docker_arch.gnu
RUN cd boxmodel_source && ./build.sh --arch docker_arch.gnu

# -----------------------------------------------------------------------------
# Stage 2: Python Dependencies Stage
# -----------------------------------------------------------------------------
FROM python:3.9-slim AS python-builder

# Install build dependencies for Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy and install Python requirements
COPY gecko_web/requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# -----------------------------------------------------------------------------
# Stage 3: Production Runtime
# -----------------------------------------------------------------------------
FROM python:3.9-slim AS production

# Labels for container registry
LABEL org.opencontainers.image.title="GECKO-A Web Interface"
LABEL org.opencontainers.image.version="3.0.6"
LABEL org.opencontainers.image.description="Atmospheric Chemistry Mechanism Generator"
LABEL org.opencontainers.image.authors="Deeksha Sharma"
LABEL org.opencontainers.image.source="https://github.com/yourusername/gecko-web"

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Fortran runtime libraries
    libgfortran5 \
    libnetcdf19 \
    libnetcdff7 \
    # Graphviz for diagram generation (system binary required)
    graphviz \
    # Health check utility
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user for security
RUN groupadd -g 1001 gecko && \
    useradd -u 1001 -g gecko -m -s /bin/bash gecko

WORKDIR /app

# Copy Python packages from builder
COPY --from=python-builder /install /usr/local

# Copy compiled Fortran binaries from builder
COPY --from=fortran-builder /build/gecko_source ./gecko_source
COPY --from=fortran-builder /build/boxmodel_source ./boxmodel_source

# Copy application code
COPY --chown=gecko:gecko gecko_web ./gecko_web
COPY --chown=gecko:gecko samples ./samples

# Configure GECKO-A paths for container environment
RUN sed -i "s|dirgecko=.*|dirgecko='/app/gecko_source/'|" gecko_source/INPUT/gecko.nml && \
    sed -i "s|dirout=.*|dirout='/app/gecko_source/RUN/OUT/'|" gecko_source/INPUT/gecko.nml

# Create data directories with correct permissions
RUN mkdir -p /data/output /data/mechanisms /data/library /data/archives \
    && chown -R gecko:gecko /data \
    && chown -R gecko:gecko /app

# Create workspace directory for Fortran (short path for 100-char limit)
RUN mkdir -p /tmp/gk && chown gecko:gecko /tmp/gk

# Set environment variables
ENV DATA_DIR=/data \
    GECKO_SOURCE_DIR=/app/gecko_source \
    BOXMODEL_SOURCE_DIR=/app/boxmodel_source \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    # Matplotlib non-interactive backend
    MPLBACKEND=Agg

# Switch to non-root user
USER gecko

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://localhost:8000/api/environment || exit 1

# Start application
CMD ["uvicorn", "gecko_web.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### 2.3 Optimized Docker Compose

Create `docker/docker-compose.optimized.yml`:

```yaml
version: '3.8'

services:
  gecko-app:
    build:
      context: ..
      dockerfile: docker/Dockerfile.optimized
      # Multi-platform builds (uncomment for CI/CD)
      # platforms:
      #   - linux/amd64
      #   - linux/arm64
    image: gecko-app:${VERSION:-latest}
    container_name: gecko-web
    ports:
      - "${PORT:-8000}:8000"
    environment:
      - DATA_DIR=/data
      - GECKO_SOURCE_DIR=/app/gecko_source
      - BOXMODEL_SOURCE_DIR=/app/boxmodel_source
    volumes:
      # Persistent job data
      - gecko-data:/data
      # Development: mount source for hot reload (comment out for production)
      # - ../gecko_web:/app/gecko_web:ro
    restart: unless-stopped
    # Resource limits
    deploy:
      resources:
        limits:
          cpus: '2.0'
          memory: 4G
        reservations:
          cpus: '0.5'
          memory: 512M
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/api/environment"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 15s
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"

volumes:
  gecko-data:
    driver: local

networks:
  default:
    name: gecko-network
```

### 2.4 Create .dockerignore

Create `docker/.dockerignore`:

```
# =============================================================================
# GECKO-A Docker Ignore File
# Reduces build context from ~1.3GB to ~65MB
# =============================================================================

# Python virtual environments (465MB+)
.venv/
venv/
__pycache__/
*.pyc
*.pyo
*.pyd
.Python
env/
.env

# Windows installer (645MB)
Gecko Windows Installation/

# Test artifacts
tests/
test_environment/
test_vis_bench/
TestGeko/
.pytest_cache/
.coverage
htmlcov/
*.egg-info/

# IDE and editor files
.idea/
.vscode/
*.swp
*.swo
*.sublime-*
.DS_Store
Thumbs.db

# Git
.git/
.gitignore
.gitattributes

# Documentation (not needed in container)
*.md
!README.md
docs/

# Local scripts (not needed in container)
start_gecko.command
compile_local.sh
rebuild_docker.sh
run_tests.sh
installers/
scripts/

# Build artifacts
*.log
logs/
tmp/
temp/
*.tmp

# Data directories (mounted as volumes)
data/

# Fortran build artifacts (rebuilt in container)
docker/gecko_source/OBJ/*.o
docker/gecko_source/OBJ/*.mod
docker/gecko_source/RUN/OUT/*
docker/boxmodel_source/OBJ/*
docker/boxmodel_source/SIMU/*/
```

### 2.5 Add Health Check Endpoint

Add to `gecko_web/main.py` (if not exists):

```python
@app.get("/health")
async def health_check():
    """Health check endpoint for container orchestration."""
    return {
        "status": "healthy",
        "version": "3.0.6",
        "gecko_available": GECKO_ROOT.exists() and (GECKO_ROOT / "RUN" / "gecko.sh").exists(),
        "boxmodel_available": BOXMODEL_ROOT.exists() and (BOXMODEL_ROOT / "prepare_simu.sh").exists(),
        "data_dir_writable": os.access(DATA_DIR, os.W_OK)
    }
```

---

## Phase 3: Cross-Platform Compatibility

### 3.1 Windows Support Files

#### 3.1.1 Create `scripts/start.bat` (Windows Batch)

```batch
@echo off
setlocal EnableDelayedExpansion

:: ============================================================================
:: GECKO-A Web Interface - Windows Launcher
:: Version 3.0.6
:: ============================================================================

echo.
echo ========================================================================
echo   GECKO-A Web Interface - Windows Launcher
echo   Version 3.0.6
echo ========================================================================
echo.

cd /d "%~dp0\.."

:: Check for Docker
where docker >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    echo [OK] Docker found
    goto :docker_mode
) else (
    echo [WARN] Docker not found - attempting local Python mode
    goto :python_mode
)

:docker_mode
echo.
echo Starting with Docker (recommended for full functionality)...
echo.

:: Check if docker-compose exists
where docker-compose >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    docker-compose -f docker/docker-compose.yml up -d
) else (
    docker compose -f docker/docker-compose.yml up -d
)

if %ERRORLEVEL% EQU 0 (
    echo.
    echo [SUCCESS] GECKO-A is running!
    echo URL: http://localhost:8000
    echo.
    echo To stop: docker-compose -f docker/docker-compose.yml down
    timeout /t 5 >nul
    start http://localhost:8000
) else (
    echo [ERROR] Docker failed to start. Trying Python mode...
    goto :python_mode
)
goto :end

:python_mode
echo.
echo Starting in LIMITED mode (UI only, no simulations)
echo For full functionality, install Docker Desktop
echo.

:: Check for Python
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Python is not installed.
    echo Please install Python 3.9+ from https://www.python.org/downloads/
    pause
    exit /b 1
)

:: Check Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYVER=%%i
echo [OK] Python %PYVER%

:: Create virtual environment if needed
if not exist "venv" (
    echo Creating virtual environment...
    python -m venv venv
)

:: Activate virtual environment
call venv\Scripts\activate.bat

:: Install dependencies
echo Installing dependencies...
pip install -q --upgrade pip
pip install -q -r gecko_web\requirements.txt

:: Set environment variables
set DATA_DIR=%cd%\data
set GECKO_SOURCE_DIR=%cd%\docker\gecko_source
set BOXMODEL_SOURCE_DIR=%cd%\docker\boxmodel_source

:: Create data directory
if not exist "data\output" mkdir data\output

:: Start server
echo.
echo [SUCCESS] Starting GECKO-A Web Interface
echo URL: http://localhost:8000
echo Press Ctrl+C to stop
echo.

start http://localhost:8000
uvicorn gecko_web.main:app --host 127.0.0.1 --port 8000

:end
endlocal
```

#### 3.1.2 Create `scripts/start.ps1` (PowerShell)

```powershell
#Requires -Version 5.1
# ============================================================================
# GECKO-A Web Interface - PowerShell Launcher
# Version 3.0.6
# ============================================================================

$ErrorActionPreference = "Stop"
$Version = "3.0.6"

function Write-Header {
    Write-Host ""
    Write-Host "========================================================================" -ForegroundColor Cyan
    Write-Host "  GECKO-A Web Interface - PowerShell Launcher" -ForegroundColor White
    Write-Host "  Version $Version" -ForegroundColor Green
    Write-Host "========================================================================" -ForegroundColor Cyan
    Write-Host ""
}

function Test-Docker {
    try {
        $null = Get-Command docker -ErrorAction Stop
        $dockerInfo = docker info 2>&1
        return $LASTEXITCODE -eq 0
    } catch {
        return $false
    }
}

function Start-DockerMode {
    Write-Host "[INFO] Starting with Docker (recommended for full functionality)" -ForegroundColor Yellow
    Write-Host ""

    Push-Location (Split-Path $PSScriptRoot -Parent)

    try {
        # Try docker compose (v2) first, then docker-compose (v1)
        $composeCmd = if (Get-Command "docker" -ErrorAction SilentlyContinue) {
            "docker compose"
        } else {
            "docker-compose"
        }

        Invoke-Expression "$composeCmd -f docker/docker-compose.yml up -d"

        if ($LASTEXITCODE -eq 0) {
            Write-Host ""
            Write-Host "[SUCCESS] GECKO-A is running!" -ForegroundColor Green
            Write-Host "URL: http://localhost:8000" -ForegroundColor Cyan
            Write-Host ""
            Write-Host "To stop: $composeCmd -f docker/docker-compose.yml down" -ForegroundColor Gray

            Start-Sleep -Seconds 3
            Start-Process "http://localhost:8000"
        } else {
            throw "Docker compose failed"
        }
    } catch {
        Write-Host "[ERROR] Docker failed. Falling back to Python mode..." -ForegroundColor Red
        Start-PythonMode
    } finally {
        Pop-Location
    }
}

function Start-PythonMode {
    Write-Host ""
    Write-Host "[WARN] Starting in LIMITED mode (UI only, no simulations)" -ForegroundColor Yellow
    Write-Host "For full functionality, install Docker Desktop" -ForegroundColor Yellow
    Write-Host ""

    Push-Location (Split-Path $PSScriptRoot -Parent)

    try {
        # Check Python
        $python = Get-Command python -ErrorAction Stop
        $pyVersion = & python --version
        Write-Host "[OK] $pyVersion" -ForegroundColor Green

        # Create venv if needed
        if (-not (Test-Path "venv")) {
            Write-Host "Creating virtual environment..." -ForegroundColor Yellow
            & python -m venv venv
        }

        # Activate venv
        & .\venv\Scripts\Activate.ps1

        # Install dependencies
        Write-Host "Installing dependencies..." -ForegroundColor Yellow
        & pip install -q --upgrade pip
        & pip install -q -r gecko_web\requirements.txt

        # Set environment
        $env:DATA_DIR = Join-Path $PWD "data"
        $env:GECKO_SOURCE_DIR = Join-Path $PWD "docker\gecko_source"
        $env:BOXMODEL_SOURCE_DIR = Join-Path $PWD "docker\boxmodel_source"

        # Create data directory
        New-Item -ItemType Directory -Force -Path "data\output" | Out-Null

        Write-Host ""
        Write-Host "[SUCCESS] Starting GECKO-A Web Interface" -ForegroundColor Green
        Write-Host "URL: http://localhost:8000" -ForegroundColor Cyan
        Write-Host "Press Ctrl+C to stop" -ForegroundColor Gray
        Write-Host ""

        Start-Process "http://localhost:8000"
        & uvicorn gecko_web.main:app --host 127.0.0.1 --port 8000

    } catch {
        Write-Host "[ERROR] Python not found. Please install Python 3.9+ from https://www.python.org" -ForegroundColor Red
        exit 1
    } finally {
        Pop-Location
    }
}

# Main execution
Write-Header

if (Test-Docker) {
    Write-Host "[OK] Docker found and running" -ForegroundColor Green
    Start-DockerMode
} else {
    Write-Host "[WARN] Docker not available" -ForegroundColor Yellow
    Start-PythonMode
}
```

#### 3.1.3 Update `scripts/start.sh` (Unix/Linux/macOS)

```bash
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
```

### 3.2 Line Endings Configuration

Create `.gitattributes`:

```gitattributes
# Auto detect text files and perform LF normalization
* text=auto eol=lf

# Force LF for shell scripts
*.sh text eol=lf
*.command text eol=lf

# Force CRLF for Windows scripts
*.bat text eol=crlf
*.cmd text eol=crlf
*.ps1 text eol=crlf

# Binary files
*.png binary
*.jpg binary
*.jpeg binary
*.gif binary
*.ico binary
*.pdf binary
*.woff binary
*.woff2 binary
*.ttf binary
*.eot binary
*.pyc binary
*.o binary
*.mod binary
*.exe binary

# Fortran source (LF for consistency)
*.f text eol=lf
*.f90 text eol=lf
*.F text eol=lf
*.F90 text eol=lf

# Data files
*.csv text eol=lf
*.json text eol=lf
*.yml text eol=lf
*.yaml text eol=lf
*.nml text eol=lf
*.dat text eol=lf
*.txt text eol=lf
```

### 3.3 Cross-Platform Path Handling

The application already uses `pathlib.Path` which handles cross-platform paths. Key patterns already in use:

```python
# In main.py - already cross-platform
BASE_DIR = Path(__file__).parent.absolute()
DATA_DIR = Path(os.getenv("DATA_DIR", "/data"))
GECKO_ROOT = Path(os.getenv("GECKO_SOURCE_DIR", "/app/gecko_source"))
```

---

## Phase 4: Dependency Management

### 4.1 Pin All Dependencies

Create `gecko_web/requirements.lock.txt`:

```txt
# GECKO-A Web Interface - Locked Dependencies
# Generated: 2026-01-28
# Python: 3.9+

# Core Web Framework
fastapi==0.109.2
uvicorn==0.27.1
starlette==0.36.3
jinja2==3.1.3
python-multipart==0.0.9

# Scientific Computing
numpy==1.26.4
scipy==1.12.0
pandas==2.2.0

# Data Formats
netCDF4==1.6.5
xarray==2024.1.1

# Visualization
matplotlib==3.8.2
Pillow==10.2.0
graphviz==0.20.1
networkx==3.2.1

# Chemistry
rdkit==2023.9.4

# PDF Generation
reportlab==4.1.0

# Testing
pytest==8.0.0
httpx==0.26.0

# Utilities
aiofiles==23.2.1
python-dotenv==1.0.1
```

### 4.2 Development Requirements

Create `gecko_web/requirements-dev.txt`:

```txt
# Development dependencies
-r requirements.txt

# Testing
pytest>=8.0.0
pytest-cov>=4.1.0
pytest-asyncio>=0.23.0
httpx>=0.26.0

# Code quality
black>=24.1.0
isort>=5.13.0
flake8>=7.0.0
mypy>=1.8.0

# Documentation
mkdocs>=1.5.0
mkdocs-material>=9.5.0
```

### 4.3 Download Wheels for Offline Installation

Create `scripts/download_wheels.sh`:

```bash
#!/bin/bash
# Download wheels for offline installation

WHEEL_DIR="./wheels"
mkdir -p "$WHEEL_DIR"

# Download for multiple platforms
pip download \
    -r gecko_web/requirements.txt \
    -d "$WHEEL_DIR" \
    --platform manylinux2014_x86_64 \
    --platform manylinux2014_aarch64 \
    --platform win_amd64 \
    --platform macosx_10_9_x86_64 \
    --platform macosx_11_0_arm64 \
    --python-version 39 \
    --only-binary=:all:

echo "Wheels downloaded to $WHEEL_DIR"
echo "Install offline: pip install --no-index --find-links=$WHEEL_DIR -r requirements.txt"
```

---

## Phase 5: GitHub Repository Setup

### 5.1 Repository Structure

```
GECKO/
├── .github/
│   ├── workflows/
│   │   ├── ci.yml
│   │   ├── docker-publish.yml
│   │   └── release.yml
│   ├── ISSUE_TEMPLATE/
│   │   ├── bug_report.md
│   │   └── feature_request.md
│   └── PULL_REQUEST_TEMPLATE.md
├── docker/
│   ├── Dockerfile
│   ├── Dockerfile.optimized
│   ├── docker-compose.yml
│   ├── docker-compose.optimized.yml
│   ├── .dockerignore
│   ├── docker_arch.gnu
│   ├── gecko_source/
│   └── boxmodel_source/
├── gecko_web/
│   ├── main.py
│   ├── requirements.txt
│   ├── requirements.lock.txt
│   └── [other modules]
├── scripts/
│   ├── start.sh
│   ├── start.bat
│   └── start.ps1
├── tests/
├── samples/
├── data/
│   └── .gitkeep
├── .dockerignore (symlink to docker/.dockerignore)
├── .gitignore
├── .gitattributes
├── .env.example
├── LICENSE
├── README.md
├── DOCUMENTATION.md
├── CHANGELOG.md
└── DEPLOYMENT_INSTRUCTIONS.md
```

### 5.2 Create .gitignore

```gitignore
# =============================================================================
# GECKO-A Web Interface - Git Ignore
# =============================================================================

# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Virtual environments
.venv/
venv/
ENV/
env/

# IDE
.idea/
.vscode/
*.swp
*.swo
*.sublime-workspace
*.sublime-project

# OS
.DS_Store
.DS_Store?
._*
.Spotlight-V100
.Trashes
ehthumbs.db
Thumbs.db

# Testing
.tox/
.nox/
.coverage
.coverage.*
htmlcov/
.pytest_cache/
nosetests.xml
coverage.xml
*.cover
.hypothesis/

# Jupyter
.ipynb_checkpoints

# Environment
.env
.env.local
.env.*.local
*.pem
*.key
secrets/

# Logs
logs/
*.log
npm-debug.log*

# Temporary files
tmp/
temp/
*.tmp
*.bak

# Data and job outputs (mounted volumes)
data/output/*
data/mechanisms/*
data/archives/*
data/library/*
!data/.gitkeep

# Fortran build artifacts
docker/gecko_source/OBJ/*.o
docker/gecko_source/OBJ/*.mod
docker/gecko_source/RUN/OUT/*
docker/boxmodel_source/OBJ/*
docker/boxmodel_source/SIMU/*/OUT/*

# Windows installer (too large for git)
Gecko Windows Installation/

# Test directories
test_environment/
test_vis_bench/
TestGeko/
```

### 5.3 Create .env.example

```env
# =============================================================================
# GECKO-A Web Interface - Environment Configuration
# Copy to .env and modify as needed
# =============================================================================

# Application port
PORT=8000

# Data storage directory
DATA_DIR=/data

# GECKO-A source directory (Docker default)
GECKO_SOURCE_DIR=/app/gecko_source

# Box Model source directory (Docker default)
BOXMODEL_SOURCE_DIR=/app/boxmodel_source

# Optional: Override workspace for Fortran jobs
# (Default: /tmp/gk for short path compatibility)
# WORKSPACE_BASE=/tmp/gk

# Optional: Job retention time in seconds (default: 24 hours)
# JOB_RETENTION_SECONDS=86400
```

### 5.4 Issue Templates

Create `.github/ISSUE_TEMPLATE/bug_report.md`:

```markdown
---
name: Bug Report
about: Report a bug in GECKO-A Web Interface
title: '[BUG] '
labels: bug
assignees: ''
---

## Bug Description
A clear and concise description of the bug.

## Environment
- **OS:** [e.g., Windows 11, macOS 14.2, Ubuntu 22.04]
- **Browser:** [e.g., Chrome 120, Firefox 121]
- **GECKO-A Version:** [e.g., 3.0.6]
- **Deployment:** [Docker / Local Python]

## Steps to Reproduce
1. Go to '...'
2. Click on '....'
3. Enter '....'
4. See error

## Expected Behavior
What you expected to happen.

## Actual Behavior
What actually happened.

## Screenshots
If applicable, add screenshots.

## Console Errors
```
Paste any browser console or server errors here
```

## Additional Context
Any other context about the problem.
```

Create `.github/ISSUE_TEMPLATE/feature_request.md`:

```markdown
---
name: Feature Request
about: Suggest a new feature for GECKO-A Web Interface
title: '[FEATURE] '
labels: enhancement
assignees: ''
---

## Feature Description
A clear and concise description of the feature.

## Use Case
Why would this feature be useful? What problem does it solve?

## Proposed Solution
How you think this could be implemented.

## Alternatives Considered
Any alternative solutions you've considered.

## Additional Context
Any other context, mockups, or references.
```

### 5.5 Pull Request Template

Create `.github/PULL_REQUEST_TEMPLATE.md`:

```markdown
## Summary
Brief description of changes.

## Type of Change
- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)
- [ ] Documentation update
- [ ] Dependency update

## Changes Made
- Change 1
- Change 2
- Change 3

## Testing
- [ ] Tests pass locally
- [ ] Docker build succeeds
- [ ] Tested on macOS
- [ ] Tested on Windows
- [ ] Tested on Linux

## Screenshots (if applicable)

## Checklist
- [ ] Code follows project style guidelines
- [ ] Self-reviewed my code
- [ ] Commented hard-to-understand areas
- [ ] Updated documentation (if needed)
- [ ] Added tests (if applicable)
- [ ] All tests pass
```

---

## Phase 6: CI/CD Workflows

### 6.1 Continuous Integration

Create `.github/workflows/ci.yml`:

```yaml
name: CI

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  # Python tests on multiple OS
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
        python-version: ['3.9', '3.10', '3.11', '3.12']
      fail-fast: false

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
          cache: 'pip'
          cache-dependency-path: gecko_web/requirements.txt

      # Install Graphviz (system dependency)
      - name: Install Graphviz (Ubuntu)
        if: runner.os == 'Linux'
        run: sudo apt-get update && sudo apt-get install -y graphviz

      - name: Install Graphviz (macOS)
        if: runner.os == 'macOS'
        run: brew install graphviz

      - name: Install Graphviz (Windows)
        if: runner.os == 'Windows'
        run: choco install graphviz -y

      - name: Install Python dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r gecko_web/requirements.txt
          pip install pytest pytest-cov pytest-asyncio

      - name: Run tests
        run: |
          pytest tests/ -v --cov=gecko_web --cov-report=xml
        env:
          DATA_DIR: ./data

      - name: Upload coverage
        if: matrix.os == 'ubuntu-latest' && matrix.python-version == '3.11'
        uses: codecov/codecov-action@v4
        with:
          file: ./coverage.xml

  # Docker build test
  docker:
    runs-on: ubuntu-latest
    needs: test

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Build Docker image
        uses: docker/build-push-action@v5
        with:
          context: .
          file: docker/Dockerfile
          push: false
          tags: gecko-app:test
          cache-from: type=gha
          cache-to: type=gha,mode=max

      - name: Test Docker container
        run: |
          docker run -d --name gecko-test -p 8000:8000 gecko-app:test
          sleep 15
          curl -f http://localhost:8000/api/environment || exit 1
          docker logs gecko-test
          docker stop gecko-test

  # Code quality
  lint:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Install linters
        run: |
          pip install flake8 black isort

      - name: Check formatting with Black
        run: black --check gecko_web/ tests/

      - name: Check imports with isort
        run: isort --check-only gecko_web/ tests/

      - name: Lint with flake8
        run: flake8 gecko_web/ tests/ --max-line-length=120 --ignore=E501,W503
```

### 6.2 Docker Publish

Create `.github/workflows/docker-publish.yml`:

```yaml
name: Docker Publish

on:
  release:
    types: [published]
  push:
    branches: [main]
    tags: ['v*']
  workflow_dispatch:

env:
  REGISTRY: ghcr.io
  IMAGE_NAME: ${{ github.repository }}

jobs:
  build-and-push:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up QEMU
        uses: docker/setup-qemu-action@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Login to GitHub Container Registry
        uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY }}/${{ env.IMAGE_NAME }}
          tags: |
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=sha,prefix=
            type=raw,value=latest,enable=${{ github.ref == 'refs/heads/main' }}

      - name: Build and push
        uses: docker/build-push-action@v5
        with:
          context: .
          file: docker/Dockerfile
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max

      - name: Generate SBOM
        uses: anchore/sbom-action@v0
        with:
          image: ${{ env.REGISTRY }}/${{ env.IMAGE_NAME }}:latest
```

### 6.3 Release Workflow

Create `.github/workflows/release.yml`:

```yaml
name: Release

on:
  push:
    tags: ['v*']

permissions:
  contents: write

jobs:
  release:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout code
        uses: actions/checkout@v4
        with:
          fetch-depth: 0

      - name: Generate changelog
        id: changelog
        uses: metcalfc/changelog-generator@v4.3.1
        with:
          myToken: ${{ secrets.GITHUB_TOKEN }}

      - name: Create Release
        uses: softprops/action-gh-release@v1
        with:
          body: |
            ## Changes
            ${{ steps.changelog.outputs.changelog }}

            ## Docker
            ```bash
            docker pull ghcr.io/${{ github.repository }}:${{ github.ref_name }}
            ```

            ## Quick Start
            ```bash
            docker run -p 8000:8000 ghcr.io/${{ github.repository }}:${{ github.ref_name }}
            ```
          draft: false
          prerelease: false
```

---

## Phase 7: Cloud Deployment

### 7.1 Railway (Recommended - Easiest)

Create `railway.toml`:

```toml
[build]
builder = "dockerfile"
dockerfilePath = "docker/Dockerfile"

[deploy]
healthcheckPath = "/api/environment"
healthcheckTimeout = 100
restartPolicyType = "on_failure"
restartPolicyMaxRetries = 3

[environments]
[environments.production]
[environments.production.variables]
PORT = "8000"
DATA_DIR = "/data"
```

### 7.2 Render

Create `render.yaml`:

```yaml
services:
  - type: web
    name: gecko-web
    env: docker
    dockerfilePath: ./docker/Dockerfile
    dockerContext: .
    healthCheckPath: /api/environment
    envVars:
      - key: PORT
        value: 8000
      - key: DATA_DIR
        value: /data
    disk:
      name: gecko-data
      mountPath: /data
      sizeGB: 10
```

### 7.3 Fly.io

Create `fly.toml`:

```toml
app = "gecko-web"
primary_region = "iad"

[build]
  dockerfile = "docker/Dockerfile"

[http_service]
  internal_port = 8000
  force_https = true
  auto_stop_machines = true
  auto_start_machines = true
  min_machines_running = 0

[http_service.concurrency]
  type = "connections"
  hard_limit = 100
  soft_limit = 80

[[http_service.checks]]
  grace_period = "30s"
  interval = "30s"
  method = "GET"
  timeout = "10s"
  path = "/api/environment"

[[vm]]
  cpu_kind = "shared"
  cpus = 2
  memory_mb = 2048

[mounts]
  source = "gecko_data"
  destination = "/data"
```

### 7.4 DigitalOcean App Platform

Create `.do/app.yaml`:

```yaml
name: gecko-web
services:
  - name: web
    dockerfile_path: docker/Dockerfile
    github:
      branch: main
      deploy_on_push: true
      repo: yourusername/gecko-web
    http_port: 8000
    instance_size_slug: professional-xs
    instance_count: 1
    routes:
      - path: /
    health_check:
      http_path: /api/environment
      initial_delay_seconds: 30
      period_seconds: 30
      timeout_seconds: 10
      success_threshold: 1
      failure_threshold: 3
    envs:
      - key: DATA_DIR
        value: /data
```

### 7.5 AWS ECS (Terraform)

Create `infrastructure/aws/main.tf`:

```hcl
# GECKO-A Web Interface - AWS ECS Deployment
# Requires: AWS CLI configured, Terraform installed

terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

variable "aws_region" {
  default = "us-east-1"
}

variable "image_tag" {
  default = "latest"
}

# ECR Repository
resource "aws_ecr_repository" "gecko" {
  name                 = "gecko-web"
  image_tag_mutability = "MUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}

# ECS Cluster
resource "aws_ecs_cluster" "gecko" {
  name = "gecko-cluster"

  setting {
    name  = "containerInsights"
    value = "enabled"
  }
}

# Task Definition
resource "aws_ecs_task_definition" "gecko" {
  family                   = "gecko-web"
  network_mode             = "awsvpc"
  requires_compatibilities = ["FARGATE"]
  cpu                      = "1024"
  memory                   = "2048"
  execution_role_arn       = aws_iam_role.ecs_execution.arn

  container_definitions = jsonencode([
    {
      name  = "gecko-web"
      image = "${aws_ecr_repository.gecko.repository_url}:${var.image_tag}"
      portMappings = [
        {
          containerPort = 8000
          hostPort      = 8000
          protocol      = "tcp"
        }
      ]
      environment = [
        { name = "DATA_DIR", value = "/data" }
      ]
      logConfiguration = {
        logDriver = "awslogs"
        options = {
          awslogs-group         = "/ecs/gecko-web"
          awslogs-region        = var.aws_region
          awslogs-stream-prefix = "ecs"
        }
      }
      healthCheck = {
        command     = ["CMD-SHELL", "curl -f http://localhost:8000/api/environment || exit 1"]
        interval    = 30
        timeout     = 10
        retries     = 3
        startPeriod = 60
      }
    }
  ])
}

# IAM Role for ECS
resource "aws_iam_role" "ecs_execution" {
  name = "gecko-ecs-execution-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "ecs-tasks.amazonaws.com"
        }
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "ecs_execution" {
  role       = aws_iam_role.ecs_execution.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

# CloudWatch Log Group
resource "aws_cloudwatch_log_group" "gecko" {
  name              = "/ecs/gecko-web"
  retention_in_days = 30
}

output "ecr_repository_url" {
  value = aws_ecr_repository.gecko.repository_url
}

output "ecs_cluster_name" {
  value = aws_ecs_cluster.gecko.name
}
```

---

## Phase 8: Execution Checklist

### Pre-Deployment Verification

#### Codebase & Configuration
- [ ] All Python imports verified working
- [ ] `requirements.txt` dependencies pinned
- [ ] `requirements.lock.txt` created
- [ ] Environment variables documented in `.env.example`
- [ ] Health check endpoint (`/health` or `/api/environment`) working

#### Docker
- [ ] Dockerfile optimized with multi-stage build
- [ ] `.dockerignore` created (reduces context from 1.3GB to ~65MB)
- [ ] `docker-compose.yml` configured
- [ ] Non-root user configured
- [ ] Health check in Dockerfile
- [ ] Local `docker build` succeeds
- [ ] Local `docker run` starts application
- [ ] API accessible at `http://localhost:8000`

#### Cross-Platform
- [ ] `start.sh` tested on macOS/Linux
- [ ] `start.bat` tested on Windows
- [ ] `start.ps1` tested on Windows PowerShell
- [ ] `.gitattributes` configured for line endings
- [ ] Path handling uses `pathlib.Path`

#### Repository
- [ ] `.gitignore` comprehensive
- [ ] `.gitattributes` configured
- [ ] `README.md` updated with Docker instructions
- [ ] `CHANGELOG.md` up to date
- [ ] Issue templates created
- [ ] PR template created

#### CI/CD
- [ ] `.github/workflows/ci.yml` created
- [ ] `.github/workflows/docker-publish.yml` created
- [ ] Tests pass in CI
- [ ] Docker builds in CI
- [ ] Multi-platform images build (amd64, arm64)

### GitHub Upload Commands

```bash
# 1. Initialize git (if not already)
cd /Users/rgpro/Desktop/GECKO
git init

# 2. Create .gitignore from template above
# 3. Create .gitattributes from template above

# 4. Add all files
git add .

# 5. Create initial commit
git commit -m "Initial commit: GECKO-A Web Interface v3.0.6

- FastAPI backend with atmospheric chemistry simulation
- Docker containerization with Fortran compilation
- Cross-platform support (macOS, Windows, Linux)
- 150+ VOC compound database
- RDKit molecular visualization
- GitHub Actions CI/CD

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"

# 6. Create GitHub repository (using gh CLI)
gh repo create gecko-web --public --description "GECKO-A Atmospheric Chemistry Web Interface"

# 7. Push to GitHub
git remote add origin https://github.com/yourusername/gecko-web.git
git branch -M main
git push -u origin main

# 8. Enable GitHub Container Registry
# Go to Settings > Packages > Enable improved container support
```

### Shareable Links After Deployment

| Resource | URL |
|----------|-----|
| **GitHub Repository** | `https://github.com/yourusername/gecko-web` |
| **Docker Image (GHCR)** | `ghcr.io/yourusername/gecko-web:latest` |
| **Railway** | `https://gecko-web.up.railway.app` |
| **Render** | `https://gecko-web.onrender.com` |
| **Fly.io** | `https://gecko-web.fly.dev` |

### Quick Start for End Users

```bash
# Option 1: Docker (recommended)
docker run -p 8000:8000 ghcr.io/yourusername/gecko-web:latest

# Option 2: Docker Compose (with persistent data)
git clone https://github.com/yourusername/gecko-web.git
cd gecko-web
docker-compose -f docker/docker-compose.yml up -d

# Option 3: Local Python (UI only, no simulations)
git clone https://github.com/yourusername/gecko-web.git
cd gecko-web
./scripts/start.sh  # macOS/Linux
# or
.\scripts\start.ps1  # Windows PowerShell
```

---

## Security Considerations

1. **No secrets in repository** - Use environment variables
2. **Non-root container user** - Security best practice
3. **Health checks enabled** - Container orchestration
4. **Image scanning** - SBOM generation in CI
5. **Dependency pinning** - Reproducible builds
6. **HTTPS enforcement** - In cloud deployments

## Performance Optimization

1. **Multi-stage Docker build** - Smaller image (~800MB vs ~2GB)
2. **Layer caching** - Faster rebuilds
3. **Alpine-based images** - Where possible
4. **Graphviz system binary** - Required for diagrams
5. **Matplotlib Agg backend** - Headless rendering
6. **Job retention policy** - 24-hour cleanup

## Monitoring Recommendations

1. **Health endpoint**: `/api/environment` or `/health`
2. **Logging**: JSON format for parsing
3. **Metrics**: Consider adding `/metrics` with Prometheus format
4. **Error tracking**: Integrate Sentry or similar

---

*Document generated: 2026-01-28*
*GECKO-A Web Interface v3.0.6*
