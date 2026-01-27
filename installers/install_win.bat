@echo off
echo === GECKO-A Installer for Windows ===

REM Check for Docker
docker --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker is not installed.
    echo Please install Docker Desktop for Windows: https://www.docker.com/products/docker-desktop/
    pause
    exit /b 1
)

REM Check if Docker daemon is running
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker daemon is not running.
    echo Please start Docker Desktop.
    pause
    exit /b 1
)

echo Building GECKO-A application...
cd /d "%~dp0.."
docker-compose -f docker/docker-compose.yml build

echo Starting GECKO-A application...
docker-compose -f docker/docker-compose.yml up -d

echo Installation complete!
echo The application is running at http://localhost:8000
start http://localhost:8000
pause
