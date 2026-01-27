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
