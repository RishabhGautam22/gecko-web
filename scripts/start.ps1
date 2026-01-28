#Requires -Version 5.1
# ============================================================================
# GECKO-A Web Interface - PowerShell Launcher
# Version 3.0.10
# ============================================================================

$ErrorActionPreference = "Stop"
$Version = "3.0.10"

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
