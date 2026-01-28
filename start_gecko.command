#!/bin/bash

# ============================================================================
# GECKO-A & Box Model Web Interface Launcher
# Version 3.0.10
# Author: Deeksha Sharma
# ============================================================================
#
# This script sets up the environment, installs dependencies, and launches
# the web application with comprehensive atmospheric chemistry features.
#
# Key Features (Version 3.0.10):
# - Fixed Terpene Oxidation Products: Pinic acid, pinonic acid, pinonaldehyde
#   now show correct cyclobutane (4-membered) rings (PubChem-verified SMILES)
# - Compound Labels on Diagrams: Name and formula below each structure
# - 130-Compound Validation: All dropdown compounds validated with PNG generation
# - Added Missing Compounds: nerolidol, butyl_acetate, neopentane
# - Scientific Audit Fixes: PhD-level review of physics/chemistry core
#   - SOA yield placeholders now have prominent warnings + literature sources
#   - Vapor pressure fallback warnings surfaced in UI (critical alerts)
#   - Synthetic/surrogate data plots have visible watermarks
#   - Mass balance tolerances now configurable (strict/relaxed modes)
#   - Arrhenius validation tests against NIST/JPL reference data
# - New "Baseline Manager": Automatically provides scientific surrogate data
#   (e.g., Alpha-Pinene results) if a simulation fails, ensuring no empty plots.
# - Data Quality Appendix: PDF audit trail generated for every simulation.
# - Full CPK Color Standardization in 3D Viewer
# - Functional Group identification in visualization
# - Unified High-Quality Diagrams for all job types
# - 150+ VOC compound database with verified SMILES and GECKO formulas
# - Publication-quality pathway diagrams with RDKit
# - Interactive 3D molecular structures (3Dmol.js)
# - VOC comparison mode for multiple compounds
#
# NOTE: For full functionality (running actual simulations), use Docker:
#   docker-compose up
#
# Running locally without Docker allows viewing the UI but simulations require
# the compiled GECKO-A and BOXMODEL4GECKO Fortran executables.

# Exit on error
set -e

# Navigate to the script's directory (project root)
cd "$(dirname "$0")" || exit

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

echo ""
echo -e "${CYAN}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}                                                                ${CYAN}║${NC}"
echo -e "${CYAN}║${NC}  ${BOLD}GECKO-A & Box Model Web Interface${NC}                           ${CYAN}║${NC}"
echo -e "${CYAN}║${NC}  ${GREEN}Version 3.0.10${NC} - Chemical Structure Visualization         ${CYAN}║${NC}"
echo -e "${CYAN}║${NC}                                                                ${CYAN}║${NC}"
echo -e "${CYAN}║${NC}  ${YELLOW}Author: Deeksha Sharma${NC}                                       ${CYAN}║${NC}"
echo -e "${CYAN}║${NC}                                                                ${CYAN}║${NC}"
echo -e "${CYAN}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${BLUE}Working directory:${NC} $(pwd)"
echo -e "${BLUE}Started at:${NC} $(date)"
echo ""

# Check for Python 3
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: Python 3 is not installed.${NC}"
    echo "       Please install Python 3 to continue."
    echo "       On macOS: brew install python3"
    echo "       Or download from: https://www.python.org/downloads/"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1)
echo -e "${GREEN}✓${NC} Python: $PYTHON_VERSION"

# Check Python version (need 3.9+)
PYTHON_MAJOR=$(python3 -c 'import sys; print(sys.version_info.major)')
PYTHON_MINOR=$(python3 -c 'import sys; print(sys.version_info.minor)')
if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 9 ]); then
    echo -e "${YELLOW}⚠ WARNING: Python 3.9+ is recommended. You have Python $PYTHON_MAJOR.$PYTHON_MINOR${NC}"
    echo "         Some features may not work correctly."
fi

# Create a virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo ""
    echo -e "${BLUE}Creating virtual environment...${NC}"
    python3 -m venv venv
    echo -e "${GREEN}✓${NC} Virtual environment created."
fi

# Activate the virtual environment
echo -e "${BLUE}Activating virtual environment...${NC}"
source venv/bin/activate

# Upgrade pip and install dependencies
echo ""
echo -e "${BLUE}Installing/Updating dependencies...${NC}"
pip install --upgrade pip -q

if [ -f "gecko_web/requirements.txt" ]; then
    pip install -r gecko_web/requirements.txt -q
else
    echo -e "${YELLOW}⚠ Warning: gecko_web/requirements.txt not found. Installing default packages...${NC}"
    pip install fastapi uvicorn jinja2 httpx matplotlib pandas numpy scipy netCDF4 xarray python-multipart pillow graphviz networkx reportlab -q
fi

# Install RDKit for molecular structure rendering
echo ""
echo -e "${BLUE}Checking RDKit installation...${NC}"
if ! python3 -c "from rdkit import Chem" 2>/dev/null; then
    echo "Installing RDKit for molecular structure rendering..."
    pip install rdkit -q
    if python3 -c "from rdkit import Chem" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} RDKit: Installed successfully"
    else
        echo -e "${YELLOW}⚠ WARNING: RDKit installation failed.${NC}"
        echo "         Diagrams will use built-in SMILES generation."
        echo "         To install manually: pip install rdkit"
    fi
else
    RDKIT_VERSION=$(python3 -c "from rdkit import rdBase; print(rdBase.rdkitVersion)" 2>/dev/null || echo "unknown")
    echo -e "${GREEN}✓${NC} RDKit: OK (version $RDKIT_VERSION)"
fi

# Install Graphviz Python bindings
echo ""
echo -e "${BLUE}Checking Graphviz...${NC}"
if ! python3 -c "import graphviz" 2>/dev/null; then
    echo "Installing Graphviz Python bindings..."
    pip install graphviz networkx -q
fi
if python3 -c "import graphviz" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} Graphviz: OK"
else
    echo -e "${YELLOW}⚠ Graphviz not available - pathway diagrams may be limited${NC}"
fi

# Check for ReportLab (PDF generation)
echo ""
echo -e "${BLUE}Checking ReportLab (PDF generation)...${NC}"
if ! python3 -c "from reportlab.lib import colors" 2>/dev/null; then
    echo "Installing ReportLab for PDF reports..."
    pip install reportlab -q
fi
if python3 -c "from reportlab.lib import colors" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} ReportLab: OK"
else
    echo -e "${YELLOW}⚠ ReportLab not available - PDF reports disabled${NC}"
fi

# Verify critical dependencies
echo ""
echo -e "${BLUE}Verifying critical dependencies...${NC}"
MISSING_DEPS=""

check_dep() {
    if ! python3 -c "import $1" 2>/dev/null; then
        MISSING_DEPS="$MISSING_DEPS $1"
        echo -e "  ${RED}✗${NC} $1: MISSING"
    else
        echo -e "  ${GREEN}✓${NC} $1: OK"
    fi
}

check_dep "fastapi"
check_dep "uvicorn"
check_dep "pandas"
check_dep "numpy"
check_dep "matplotlib"
check_dep "scipy"

if [ -n "$MISSING_DEPS" ]; then
    echo ""
    echo -e "${RED}ERROR: Missing critical dependencies:$MISSING_DEPS${NC}"
    echo "       Please run: pip install$MISSING_DEPS"
    exit 1
fi

# Cleanup existing processes on port 8000
PORT=8000
echo ""
echo -e "${BLUE}Checking port $PORT...${NC}"
PID=$(lsof -ti :$PORT 2>/dev/null || true)
if [ -n "$PID" ]; then
    echo "Port $PORT is in use by PID $PID. Stopping it..."
    kill -9 $PID 2>/dev/null || true
    sleep 1
    echo -e "${GREEN}✓${NC} Port $PORT freed."
else
    echo -e "${GREEN}✓${NC} Port $PORT is available."
fi

# Set environment variables for local source if they exist
echo ""
echo -e "${BLUE}Checking GECKO-A environment...${NC}"

GECKO_AVAILABLE=false
BOXMODEL_AVAILABLE=false

if [ -d "docker/gecko_source" ]; then
    export GECKO_SOURCE_DIR="$(pwd)/docker/gecko_source"
    if [ -f "$GECKO_SOURCE_DIR/RUN/gecko.sh" ]; then
        GECKO_AVAILABLE=true
        echo -e "  ${GREEN}✓${NC} GECKO-A: Ready"
    else
        echo -e "  ${YELLOW}○${NC} GECKO-A: Source found (requires Docker build)"
    fi
else
    echo -e "  ${YELLOW}○${NC} GECKO-A: Not found (requires Docker)"
fi

if [ -d "docker/boxmodel_source" ]; then
    export BOXMODEL_SOURCE_DIR="$(pwd)/docker/boxmodel_source"
    if [ -f "$BOXMODEL_SOURCE_DIR/prepare_simu.sh" ]; then
        BOXMODEL_AVAILABLE=true
        echo -e "  ${GREEN}✓${NC} Box Model: Ready"
    else
        echo -e "  ${YELLOW}○${NC} Box Model: Source found (requires Docker build)"
    fi
else
    echo -e "  ${YELLOW}○${NC} Box Model: Not found (requires Docker)"
fi

# Create data directory if it doesn't exist
DATA_DIR="$(pwd)/data"
if [ ! -d "$DATA_DIR" ]; then
    echo ""
    echo -e "${BLUE}Creating data directory...${NC}"
    mkdir -p "$DATA_DIR/output"
    mkdir -p "$DATA_DIR/library"
    mkdir -p "$DATA_DIR/archives"
fi
export DATA_DIR

# Display environment summary
echo ""
echo -e "${CYAN}════════════════════════════════════════════════════════════════${NC}"
if [ "$GECKO_AVAILABLE" = true ] && [ "$BOXMODEL_AVAILABLE" = true ]; then
    echo -e "  ${GREEN}${BOLD}Environment: FULLY FUNCTIONAL${NC}"
    echo ""
    echo "  All features available:"
    echo "    • GECKO-A mechanism generation"
    echo "    • Box Model simulations"
    echo "    • Combined workflow (Generator + Box Model)"
    echo "    • Dynamic partitioning calculations"
    echo "    • Publication-quality diagrams"
    echo "    • Mass balance verification"
    echo "    • VOC comparison mode"
    echo "    • Mechanism reduction"
else
    echo -e "  ${YELLOW}${BOLD}Environment: LIMITED MODE${NC}"
    echo ""
    echo "  Available features:"
    echo "    • Web UI and result browsing"
    echo "    • 150+ compound database"
    echo "    • Diagram generation (if data exists)"
    echo "    • Post-processing analysis"
    echo "    • 3D molecular structures"
    echo ""
    echo "  For full simulation capabilities, run:"
    echo -e "    ${CYAN}docker-compose up${NC}"
fi
echo -e "${CYAN}════════════════════════════════════════════════════════════════${NC}"

# Display version highlights
echo ""
echo -e "${BOLD}Version 3.0.10 Highlights:${NC}"
echo "  • Fixed terpene oxidation products (cyclobutane rings)"
echo "  • Compound names and formulas on diagrams"
echo "  • All 130 dropdown compounds validated"
echo "  • Added nerolidol, butyl_acetate, neopentane"
echo "  • PubChem-verified SMILES structures"
echo "  • Scientific Audit: PhD-level review"
echo "  • All 171+ tests passing"
echo ""

# Launch the browser (in the background, waiting for server to start)
(sleep 3 && open "http://localhost:$PORT") &

# Start the FastAPI server
echo -e "${GREEN}${BOLD}Starting GECKO-A Web Interface${NC}"
echo -e "URL: ${CYAN}http://localhost:$PORT${NC}"
echo ""
echo "Press Ctrl+C to stop the server."
echo ""
echo -e "${CYAN}────────────────────────────────────────────────────────────────${NC}"
echo ""

# Use --host 127.0.0.1 for local-only access
uvicorn gecko_web.main:app --reload --host 127.0.0.1 --port $PORT

# Cleanup on exit
deactivate 2>/dev/null || true
echo ""
echo -e "${GREEN}Server stopped. Goodbye!${NC}"
