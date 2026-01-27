@echo off
echo Starting GECKO-A...
cd /d "%~dp0.."

REM Check if Docker daemon is running
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker daemon is not running.
    echo Please start Docker Desktop.
    pause
    exit /b 1
)

docker-compose -f docker/docker-compose.yml up -d
echo Opening browser...
start http://localhost:8000
