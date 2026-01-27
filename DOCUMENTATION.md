# GECKO-A Web Interface - Complete Documentation

**Author: Deeksha Sharma**

This document combines all project documentation into a single reference.

---

# Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Project Structure](#project-structure)
4. [Installation](#installation)
5. [API Reference](#api-reference)
6. [Architecture](#architecture)
7. [Chemical Database](#chemical-database)
8. [Advanced Workflows](#advanced-workflows)
9. [Scientific References](#scientific-references)
10. [Troubleshooting](#troubleshooting)
11. [GECKO-A Source Documentation](#gecko-a-source-documentation)
12. [Box Model Source Documentation](#box-model-source-documentation)

---

# Overview

This project provides a web-based interface for running GECKO-A and Box Model simulations. It is packaged using Docker to ensure consistency across macOS and Windows.

## Current Version: 3.0.5 - High-Quality Visuals & Fixes

### v3.0.5 Changes
- **Fixed RDKit Visualization**: Resolved font handling bugs that caused diagrams to revert to text boxes.
- **Best-in-Class Layouts**: Enabled RDKit's CoordGen engine for superior 2D molecular layouts.
- **Enhanced Readability**: Increased bond thickness and ensured stereo-annotation is visible.
- **Fixed Isopentane BUG**: Modified SMILES parser to prioritize specific formulas.

## Version 3.0.4 - Simplified & Stabilized

### v3.0.4 Changes
- **Removed Combined Workflow**: The feature was unreliable and has been removed
- **Fixed branched alkane display**: Isopentane, neopentane, isobutane now show correct structures
- **Fixed pathway diagram chaos**: Reduced node count and improved graphviz layout settings

### v3.0.3 Fixes
- Simplified Combined Workflow (now removed in 3.0.4)
- Fixed GECKO formulas for cyclic compounds
- Improved 3D structure generation for radical species

### v3.0.2 Fixes
- Fixed Combined Workflow failures (mechanism directory passing)
- Fixed Box Model failures for many compounds
- Fixed 3D Model positioning
- Fixed 3D Model loading failures

### v3.0.1 Fixes
- Fixed Combined Workflow pipeline
- Fixed 3D Structure Viewer
- Fixed Diagram Regeneration with graphviz
- Added interactive UI components for all advanced features

## Version 3.0.0 - Complete Scientific Overhaul

This major release implements 15 comprehensive improvements:

1. **Expanded VOC Database** - 150+ compounds with verified SMILES
2. **Mass Balance Verification** - Automatic atom conservation checking
3. **Combined Workflow** - Unified Generator + Box Model execution
4. **Publication-Quality Diagrams** - RDKit-rendered molecular structures
5. **Reaction Type Color Coding** - Visual distinction for different reaction types
6. **Branching Ratios** - Displayed on pathway diagram edges
7. **Variable Legend** - X notation for common functional groups
8. **Extended GECKO-A Options** - Full control over generator parameters
9. **Multiple Layout Algorithms** - Hierarchical, radial, force-directed
10. **PDF Report Generation** - Comprehensive reports with embedded plots
11. **Interactive 3D Structures** - 3Dmol.js visualization
12. **VOC Comparison Mode** - Compare multiple VOCs side-by-side
13. **Mechanism Reduction Tools** - DRGEP, lumping, PFA methods
14. **Reaction Kinetics Database** - IUPAC/JPL rate constants
15. **Stereochemistry Handling** - E/Z isomer detection

---

# Features

## Core Features
- **Real GECKO-A Integration**: Generates chemical mechanisms using the official GECKO-A Fortran code
- **Real Box Model Integration**: Simulates time-series data using BOXMODEL4GECKO
- **150+ Compound Database**: Verified SMILES, molecular properties, and reaction kinetics
- **Interactive Reaction Trees**: Cytoscape.js visualization with molecular structures
- **Publication-Quality Diagrams**: RDKit + Graphviz hierarchical pathway diagrams

## New in v3.0.0
- **3D Molecular Structures**: Interactive 3Dmol.js visualization
- **Multiple Export Formats**: KPP, MCM, and FACSIMILE mechanism formats
- **Scientific Plots**: VBS distributions, SOA yield curves, partitioning summaries
- **Mass Balance Verification**: Automatic atom conservation checking
- **VOC Comparison**: Side-by-side comparison of multiple compounds
- **Mechanism Reduction**: DRGEP, lumping, and sensitivity-based reduction
- **PDF Reports**: Comprehensive reports with embedded plots
- **Web Interface**: Easy-to-use UI for submitting jobs and viewing results
- **Dockerized**: Runs consistently on any platform supporting Docker

---

# Project Structure

```
GECKO/
├── gecko_web/                   # Python web application (FastAPI)
│   ├── main.py                  # Backend API with extended endpoints
│   ├── mechanism_diagram.py     # Mechanism parsing & diagram generation
│   ├── pathway_visualizer.py    # RDKit + Graphviz publication-quality diagrams
│   ├── enhanced_visualizer.py   # Advanced visualization with reaction colors
│   ├── postprocessing.py        # Dynamic partitioning & Fortran parser
│   ├── reaction_tree.py         # Reaction pathway analysis & SMILES conversion
│   ├── mass_balance.py          # Atom conservation verification
│   ├── combined_workflow.py     # Unified Generator + Box Model workflow
│   ├── chemdata/                # Chemical database package
│   │   ├── __init__.py          # Package exports
│   │   ├── compound_database.py # 150+ verified compounds
│   │   ├── voc_categories.py    # VOC categorization
│   │   └── reaction_data.py     # Reaction kinetics database
│   ├── templates/               # HTML templates
│   └── static/                  # CSS, JavaScript assets
│       ├── css/
│       │   └── styles.css       # Main stylesheet
│       └── js/
│           ├── api.js           # API client
│           ├── structure3d.js   # 3D visualization module
│           ├── reaction-tree.js # Cytoscape visualization
│           └── app.js           # Main application logic
├── docker/                      # Docker configuration
│   ├── Dockerfile               # Container definition
│   ├── docker-compose.yml       # Container orchestration
│   ├── gecko_source/            # GECKO-A Fortran source
│   └── boxmodel_source/         # BOXMODEL4GECKO source
├── tests/                       # Unit tests
├── data/                        # Simulation results (mounted volume)
├── start_gecko.command          # Quick-start launcher (macOS)
├── README.md                    # Quick start guide
├── CHANGELOG.md                 # Version history
└── DOCUMENTATION.md             # This file
```

---

# Installation

## Requirements

- Docker Desktop must be installed and running
- GECKO-A source code in `docker/gecko_source`
- BOXMODEL4GECKO source code in `docker/boxmodel_source`
- Python 3.9+ (for local development)

## Dependencies

### Python Packages
```bash
pip install fastapi uvicorn jinja2 httpx matplotlib pandas numpy scipy \
    netCDF4 xarray python-multipart rdkit graphviz networkx reportlab pillow
```

### System Requirements
```bash
# macOS - Install graphviz for diagram generation
brew install graphviz

# Linux
apt-get install graphviz
```

## Quick Start (macOS)

Double-click `start_gecko.command` in the project root. This will:
1. Create a Python virtual environment
2. Install all dependencies (including RDKit, Graphviz, ReportLab)
3. Launch the web interface at http://localhost:8000

## Docker (Full Functionality)

For complete simulation capabilities:

```bash
docker-compose up
```

The interface will be available at http://localhost:8000

## Local Development

1. Create virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

2. Install dependencies:
   ```bash
   pip install -r gecko_web/requirements.txt
   ```

3. Start the server:
   ```bash
   uvicorn gecko_web.main:app --reload
   ```

4. Access the interface at http://localhost:8000

---

# API Reference

## Core Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Web interface |
| `/api/environment` | GET | Check GECKO-A/BoxModel availability |
| `/api/jobs` | POST | Create new job |
| `/api/jobs` | GET | List all jobs |
| `/api/jobs/{id}` | GET | Get job status |
| `/api/jobs/{id}/results` | GET | Get job results |
| `/api/jobs/{id}/archive` | POST | Archive completed job |
| `/api/jobs/{id}` | DELETE | Delete job |

## Compound Database Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/compounds` | GET | List all compounds with properties |
| `/api/compounds/categories` | GET | Get VOC category hierarchy |
| `/api/compounds/{name}` | GET | Get detailed compound info |
| `/api/compounds/search/{query}` | GET | Search compounds |
| `/api/compounds/category/{cat}` | GET | Get compounds by category |

## Reaction Kinetics Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/kinetics/{name}` | GET | Get rate constants at temperature |
| `/api/kinetics/{name}/lifetime` | GET | Calculate atmospheric lifetime |

## Advanced Workflow Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/workflow/combined` | POST | Combined Generator + Box Model |
| `/api/workflow/comparison` | POST | Compare multiple VOCs |
| `/api/mechanism/reduce` | POST | Reduce mechanism size |
| `/api/jobs/{id}/mass-balance` | GET | Verify mass balance |
| `/api/jobs/{id}/visualize` | POST | Generate enhanced diagram |

## 3D Structure Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/structure/{name}/3d` | GET | Get 3D structure (MOL/SDF/XYZ/PDB) |
| `/api/jobs/{id}/structures/3d` | GET | Get all structures for job |

## Creating Jobs

### Basic Job
```bash
curl -X POST http://localhost:8000/api/jobs \
  -H "Content-Type: application/json" \
  -d '{
    "voc_name": "alpha_pinene",
    "job_type": "boxmodel",
    "scenario": {
      "temperature_k": 298.0,
      "max_generations": 2,
      "initial_o3_ppb": 40.0,
      "seed_aerosol_ug_m3": 10.0
    }
  }'
```

### Combined Workflow
```bash
curl -X POST http://localhost:8000/api/workflow/combined \
  -H "Content-Type: application/json" \
  -d '{
    "voc_name": "limonene",
    "run_generator": true,
    "run_boxmodel": true,
    "generate_diagrams": true,
    "generate_pdf_report": true,
    "verify_mass_balance": true,
    "generator_options": {
      "vapor_pressure_method": "nannoolal",
      "max_generations": 3,
      "enable_autoxidation": false
    }
  }'
```

### VOC Comparison
```bash
curl -X POST http://localhost:8000/api/workflow/comparison \
  -H "Content-Type: application/json" \
  -d '{
    "voc_names": ["alpha_pinene", "limonene", "isoprene"],
    "compare_mechanism_size": true,
    "compare_soa_yield": true
  }'
```

---

# Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                        FastAPI Application                           │
├─────────────────────────────────────────────────────────────────────┤
│  JobManager            │  WorkspaceManager                          │
│  - Thread-safe state   │  - Isolated workspaces per job             │
│  - JSON persistence    │  - Copies GECKO executables                │
│  - UUID job IDs        │  - Updates gecko.nml paths                 │
├─────────────────────────────────────────────────────────────────────┤
│  chemdata Package      │  Visualization Modules                     │
│  - 150+ compounds      │  - pathway_visualizer.py                   │
│  - VOC categories      │  - enhanced_visualizer.py                  │
│  - Reaction kinetics   │  - mechanism_diagram.py                    │
├─────────────────────────────────────────────────────────────────────┤
│  Processing Modules    │  Workflow Modules                          │
│  - postprocessing.py   │  - combined_workflow.py                    │
│  - mass_balance.py     │  - Mechanism reduction                     │
│  - reaction_tree.py    │  - VOC comparison                          │
├─────────────────────────────────────────────────────────────────────┤
│  Output Generators                                                   │
│  - Cytoscape JSON  - KPP format    - Static diagrams  - PDF reports │
│  - MCM format      - FACSIMILE     - VBS plots        - 3D structs  │
└─────────────────────────────────────────────────────────────────────┘
```

---

# Chemical Database

## Compound Categories

The database contains 150+ compounds organized into categories:

### Alkanes
- Linear (methane to hexadecane)
- Branched (isobutane, neopentane, etc.)
- Cyclic (cyclopentane, cyclohexane, methylcyclohexane)

### Alkenes
- Simple (ethene, propene, 1-butene)
- Dienes (isoprene, 1,3-butadiene)
- Cyclic (cyclopentene, cyclohexene)

### Aromatics
- BTEX (benzene, toluene, ethylbenzene, xylenes)
- PAHs (naphthalene, etc.)
- Phenols

### Terpenes
- Monoterpenes (alpha-pinene, beta-pinene, limonene, myrcene, etc.)
- Sesquiterpenes (beta-caryophyllene, etc.)
- Oxygenated terpenes (linalool, geraniol)

### Oxygenated VOCs
- Alcohols (methanol to butanol)
- Aldehydes (formaldehyde, acetaldehyde)
- Ketones (acetone, MEK)
- Carboxylic acids

## Data Fields

Each compound includes:
- **SMILES**: Validated molecular structure
- **GECKO Formula**: Native GECKO-A notation
- **Molecular Weight**: g/mol
- **Vapor Pressure**: Pa at 298K
- **Rate Constants**: OH, O3, NO3 at 298K
- **CAS Number**: Chemical identifier
- **InChI**: IUPAC identifier

---

# Advanced Workflows

## Combined Workflow

Runs the complete GECKO-A pipeline:
1. Mechanism Generation
2. Box Model Simulation
3. Post-processing Analysis
4. Diagram Generation
5. Mass Balance Verification
6. PDF Report Generation

## VOC Comparison

Compare multiple VOCs on:
- Mechanism size (species, reactions)
- Product distribution (carbon numbers)
- SOA yields (estimated)
- Radical budget

## Mechanism Reduction

Four reduction methods available:

### DRGEP (Directed Relation Graph with Error Propagation)
Based on Pepiot-Desjardins & Pitsch (2008). Uses path-based importance analysis to identify essential species.

### Lumping
Groups similar species by carbon number and functional groups. Useful for reducing complexity while preserving key chemistry.

### PFA (Path Flux Analysis)
Based on Sun et al. (2010). Uses flux-weighted species importance.

### Sensitivity-based
Simple connectivity analysis for quick reduction.

---

# Scientific References

- **Pankow, J.F. (1994)**: An absorption model of gas/particle partitioning of organic compounds in the atmosphere. *Atmos. Environ.* 28, 185-188.
- **Odum, J.R. et al. (1996)**: Gas/particle partitioning and secondary organic aerosol yields. *Environ. Sci. Technol.* 30, 2580-2585.
- **Donahue, N.M. et al. (2006)**: Coupled partitioning, dilution, and chemical aging of semivolatile organics. *Environ. Sci. Technol.* 40, 2635-2643.
- **Pepiot-Desjardins & Pitsch (2008)**: An efficient error-propagation-based reduction method. *Combust. Flame* 154, 67-81.
- **Sun et al. (2010)**: Path flux analysis for reduction of detailed mechanisms. *Combust. Flame* 157, 1298-1307.
- **Nannoolal et al. (2008)**: Vapor pressure estimation methods. *Fluid Phase Equilibria* 269, 117-133.

---

# Troubleshooting

## Port 8000 Already in Use
The start script automatically frees port 8000, but if issues persist:
```bash
lsof -ti :8000 | xargs kill -9
```

## RDKit Installation Issues
If RDKit fails to install via pip, try conda:
```bash
conda install -c conda-forge rdkit
```

## Docker Build Fails
Ensure Docker Desktop is running and has sufficient resources allocated (4GB+ RAM recommended).

## Jobs Stuck in "Running" State
Check the job logs via `/api/jobs/{id}` endpoint. Common issues:
- GECKO-A executable not found
- Insufficient disk space
- Invalid VOC name

## 3D Structure Generation Fails
Ensure RDKit is properly installed:
```python
from rdkit import Chem
from rdkit.Chem import AllChem
print("RDKit OK")
```

## Diagram Regeneration Fails
Both system graphviz and Python bindings are required:
```bash
# Install system graphviz
brew install graphviz      # macOS
apt-get install graphviz   # Linux

# Install Python bindings
pip install graphviz
```

Verify installation:
```python
import graphviz
print("Graphviz OK")
```

## Mass Balance Errors
Small imbalances (<1%) are normal due to:
- Floating point precision
- Incomplete mechanism parsing
- Missing inorganic species

## Viewing Downloaded 3D Structure Files
For MOL, SDF, XYZ files, use molecular visualization software:
- **Avogadro** (free, recommended): https://avogadro.cc/
- **PyMOL**: https://pymol.org/
- **UCSF ChimeraX**: https://www.rbvi.ucsf.edu/chimera/
- **Online**: MolView.org (paste file content)

---

# GECKO-A Source Documentation

A detailed wiki is available at: https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/gecko-a/-/wikis/home

## Repository Structure

### DATA
Contains all input data such as experimental rate constants, inorganic chemical schemes, forced reactions, etc.

### INPUT
Contains input files: `cheminput.dat` (precursors) and `gecko.nml` (namelist parameters).

### LIB
Contains all GECKO-A Fortran 90 source code.

### OBJ
Compilation directory. Use `make clean` then `make` to compile.

### RUN
Run directory. Execute simulations using `gecko.sh`. Results are in `OUT` folder.

## GECKO-A Credits

### Current Developers
- Bernard Aumont (UPEC/LISA)
- Marie Camredon (UPEC/LISA)
- Julia Lee-Taylor (NCAR/ACOM)
- Richard Valorso (CNRS/LISA)

---

# Box Model Source Documentation

A detailed wiki is available at: https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/boxmodel4gecko/-/wikis/home

## Repository Structure

### CHEMDAT
Contains chemical mechanisms, saturation vapor pressure files, etc.

### INPUT
Contains input files such as photolysis tables.

### INTERP
Contains the interpreter for translating chemical mechanisms to binary.

### LIBSRC
Contains all Fortran subroutines.

### OBJ
Compilation directory.

### SIMU
Simulation folder.

### PROG
Contains main boxmodel program.

## Compilation
Requirements: Fortran compiler (e.g., gfortran), netCDF library.

Use `build.sh` script with `--arch` option specifying your architecture file.

## Preparing a Simulation
```bash
./prepare_simu.sh --import_mech nameforyourmechanism --geckooutdir PATH/TO/GECKO/RUN/OUT
```

## Running a Simulation
Go to the simulation folder and run `start.sh`.

## Box Model Credits

### Current Developers
- Bernard Aumont (UPEC/LISA)
- Marie Camredon (UPEC/LISA)
- Julia Lee-Taylor (NCAR/ACOM)
- Richard Valorso (CNRS/LISA)

---

*Documentation for GECKO-A Web Interface v3.0.1*
*Author: Deeksha Sharma*
