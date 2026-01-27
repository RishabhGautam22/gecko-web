# GECKO-A & Box Model Interface

**Author: Deeksha Sharma**

This project provides a web-based interface for running GECKO-A and Box Model simulations. It is packaged using Docker to ensure consistency across macOS and Windows.

---

## 🚀 Quick Start

### Run on Cloud (Fly.io) - Recommended
```bash
# One-time setup
fly auth login
fly launch --no-deploy
fly volumes create gecko_data --region iad --size 10
fly deploy

# Access at: https://your-app-name.fly.dev
```
See [FLY_DEPLOYMENT.md](FLY_DEPLOYMENT.md) for detailed instructions.

### Run Locally (macOS)
```bash
# Double-click or run:
./start_gecko.command
# Opens http://localhost:8000
```

### Run with Docker (Full Functionality)
```bash
docker-compose -f docker/docker-compose.yml up
# Opens http://localhost:8000
```

---

## 📖 Documentation

| Document | Purpose |
|----------|---------|
| [FLY_DEPLOYMENT.md](FLY_DEPLOYMENT.md) | Deploy to Fly.io (cloud) |
| [DEPLOYMENT_INSTRUCTIONS.md](DEPLOYMENT_INSTRUCTIONS.md) | Complete deployment guide |
| [CLAUDE.md](CLAUDE.md) | Instructions for VS Code / Claude Code |
| [DOCUMENTATION.md](DOCUMENTATION.md) | Full application documentation |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

---

## Version 3.0.6 - Visualization & Workflow Consistency (Latest)

### Recent Changes (v3.0.6)
- **Unified Diagrams**: Both "Generator" and "Box Model" jobs now produce high-quality, full-color oxidation pathway diagrams.
- **Enhanced 3D Viewer**: Now uses standard CPK colors (Dark Carbon) and displays functional group names (e.g., "Nitrate", "PAN") alongside codes.
- **Fixed Workflow Bugs**: Resolved issues with invisible diagrams (`.png.png` double extension) and broken download links.

## Version 3.0.5 - Visual Perfection

### Recent Changes (v3.0.5)
- **Fixed RDKit Visualization**: Resolved font handling bugs that caused diagrams to revert to "boxes".
- **Best-in-Class Layouts**: Enabled RDKit's CoordGen engine for superior 2D molecular layouts.
- **Enhanced Readability**: Increased bond thickness (2.5x) and ensured stereo-annotation (E/Z, R/S) is visible.
- **Fixed Isopentane BUG**: Modified SMILES parser to prioritize specific formulas over generic GECKO codes.

## Version 3.0.4 - Simplified & Stabilized

### Recent Changes (v3.0.4)
- **Removed Combined Workflow**: Feature was unreliable and has been removed entirely
- **Fixed branched alkane display**: Isopentane, neopentane, etc. now show correct structures
- **Fixed pathway diagrams**: Reduced node count and improved graphviz layout for cleaner diagrams
- **Simplified UI**: Job type now only shows Generator, Box Model, and Comparison modes

### Previous Fixes (v3.0.3)
- Fixed GECKO formulas for cyclic compounds
- Improved 3D structure generation for radical species

### Previous Fixes (v3.0.2)
- Fixed Combined Workflow failures (mechanism directory passing)
- Fixed Box Model failures for many compounds
- Fixed 3D Model positioning
- Fixed 3D Model loading failures for radical species

### Previous Fixes (v3.0.1)
- Fixed Combined Workflow pipeline (was failing due to missing config attributes)
- Fixed 3D Structure Viewer (now properly loads and displays molecules)
- Fixed Diagram Regeneration (graphviz integration corrected)
- Added interactive UI for VOC Comparison, Mechanism Reduction, and Layout selection
- Added Mass Balance verification display

## Version 3.0.0 - Complete Scientific Overhaul

This major release implements 15 comprehensive improvements to scientific accuracy, visualization, and usability:

### New Features

1. **Expanded VOC Database** - 150+ compounds with verified SMILES, molecular properties, and reaction rate constants
2. **Mass Balance Verification** - Automatic atom conservation checking (C, H, O, N, S) for generated mechanisms
3. **Combined Workflow** - Unified Generator + Box Model execution with structured output directories
4. **Publication-Quality Diagrams** - RDKit-rendered molecular structures with Graphviz layouts
5. **Reaction Type Color Coding** - Visual distinction for OH addition (blue), O3 ozonolysis (orange), NO3 (purple), photolysis (yellow)
6. **Branching Ratios** - Displayed on pathway diagram edges
7. **Variable Legend** - X notation for common functional groups (-OOH, -ONO2, -OH)
8. **Extended GECKO-A Options** - Vapor pressure methods, reaction channels, stereochemistry handling
9. **Multiple Layout Algorithms** - Hierarchical, radial, and force-directed layouts
10. **PDF Report Generation** - Comprehensive reports with embedded plots
11. **Interactive 3D Structures** - 3Dmol.js integration for molecular visualization
12. **VOC Comparison Mode** - Compare multiple VOCs side-by-side
13. **Mechanism Reduction Tools** - DRGEP, lumping, PFA, and sensitivity-based reduction
14. **Reaction Kinetics Database** - Arrhenius and Troe parameters from IUPAC/JPL
15. **Stereochemistry Handling** - E/Z isomer detection and display

See CHANGELOG.md for full details.

## Project Structure

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
│   │   ├── compound_database.py # 150+ verified compounds
│   │   ├── voc_categories.py    # VOC categorization
│   │   └── reaction_data.py     # Reaction kinetics database
│   ├── templates/               # HTML templates
│   └── static/                  # CSS, JavaScript assets
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
└── start_gecko.command          # Quick-start launcher (macOS)
```

## Features

- **Real GECKO-A Integration**: Generates chemical mechanisms using the official GECKO-A Fortran code
- **Real Box Model Integration**: Simulates time-series data using BOXMODEL4GECKO
- **150+ Compound Database**: Verified SMILES, molecular properties, and reaction kinetics
- **Interactive Reaction Trees**: Cytoscape.js visualization with molecular structures
- **Publication-Quality Diagrams**: RDKit + Graphviz hierarchical pathway diagrams
- **3D Molecular Structures**: Interactive 3Dmol.js visualization
- **Multiple Export Formats**: KPP, MCM, and FACSIMILE mechanism formats
- **Scientific Plots**: VBS distributions, SOA yield curves, partitioning summaries
- **Mass Balance Verification**: Automatic atom conservation checking
- **VOC Comparison**: Side-by-side comparison of multiple compounds
- **Mechanism Reduction**: DRGEP, lumping, and sensitivity-based reduction
- **PDF Reports**: Comprehensive reports with embedded plots
- **Web Interface**: Easy-to-use UI for submitting jobs and viewing results
- **Dockerized**: Runs consistently on any platform supporting Docker

## Requirements

- Docker Desktop must be installed and running
- GECKO-A source code in `docker/gecko_source`
- BOXMODEL4GECKO source code in `docker/boxmodel_source`
- Python 3.9+ (for local development)

### Dependencies

- **RDKit**: Molecular structure rendering, SMILES validation, and 3D structure generation
- **Graphviz**: Publication-quality pathway diagrams (requires both system package and Python bindings)
- **ReportLab**: PDF report generation
- **3Dmol.js**: Interactive 3D structures (loaded from CDN)

```bash
# Install Python packages
pip install rdkit graphviz networkx reportlab

# Install system graphviz (required for diagram generation)
brew install graphviz      # macOS
# apt-get install graphviz # Linux
```

## How to Install & Run

### Quick Start (macOS)

Double-click `start_gecko.command` in the project root. This will:
1. Create a Python virtual environment
2. Install all dependencies (including RDKit, Graphviz, ReportLab)
3. Launch the web interface at http://localhost:8000

### Docker (Full Functionality)

For complete simulation capabilities:

```bash
docker-compose up
```

The interface will be available at http://localhost:8000

### Local Development

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

## API Endpoints

### Core Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Web interface |
| `/api/environment` | GET | Check GECKO-A/BoxModel availability |
| `/api/jobs` | POST | Create new job |
| `/api/jobs` | GET | List all jobs |
| `/api/jobs/{id}` | GET | Get job status |
| `/api/jobs/{id}/results` | GET | Get job results |

### Compound Database

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/compounds` | GET | List all compounds |
| `/api/compounds/categories` | GET | Get VOC categories |
| `/api/compounds/{name}` | GET | Get compound details |
| `/api/compounds/search/{query}` | GET | Search compounds |

### Reaction Kinetics

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/kinetics/{name}` | GET | Get rate constants |
| `/api/kinetics/{name}/lifetime` | GET | Calculate atmospheric lifetime |

### Advanced Workflows

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/workflow/combined` | POST | Combined Generator + Box Model |
| `/api/workflow/comparison` | POST | Compare multiple VOCs |
| `/api/mechanism/reduce` | POST | Reduce mechanism size |
| `/api/jobs/{id}/mass-balance` | GET | Verify mass balance |
| `/api/jobs/{id}/visualize` | POST | Generate enhanced diagram |

### 3D Structures

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/structure/{name}/3d` | GET | Get 3D molecular structure |
| `/api/jobs/{id}/structures/3d` | GET | Get all 3D structures for job |

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
    "verify_mass_balance": true
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

## Scientific References

- **Pankow, J.F. (1994)**: An absorption model of gas/particle partitioning of organic compounds in the atmosphere. *Atmos. Environ.* 28, 185-188.
- **Odum, J.R. et al. (1996)**: Gas/particle partitioning and secondary organic aerosol yields. *Environ. Sci. Technol.* 30, 2580-2585.
- **Donahue, N.M. et al. (2006)**: Coupled partitioning, dilution, and chemical aging of semivolatile organics. *Environ. Sci. Technol.* 40, 2635-2643.
- **Pepiot-Desjardins & Pitsch (2008)**: DRGEP mechanism reduction method. *Combust. Flame* 154, 67-81.

## Troubleshooting

### Port 8000 Already in Use
The start script automatically frees port 8000, but if issues persist:
```bash
lsof -ti :8000 | xargs kill -9
```

### RDKit Installation Issues
If RDKit fails to install via pip, try conda:
```bash
conda install -c conda-forge rdkit
```

### 3D Structure Viewer Shows Black/Empty
- Ensure RDKit is installed: `pip install rdkit`
- Check browser console for errors
- Try selecting a different species from the dropdown

### Diagram Regeneration Fails
- Install system graphviz: `brew install graphviz` (macOS) or `apt-get install graphviz` (Linux)
- Install Python graphviz: `pip install graphviz`
- Restart the server after installation

### Combined Workflow Fails
- Check job logs for specific stage failures
- Ensure GECKO-A source is available (Docker mode recommended)
- Verify all required directories exist in `data/output/`

### Docker Build Fails
Ensure Docker Desktop is running and has sufficient resources allocated (4GB+ RAM recommended).

### Jobs Stuck in "Running" State
Check the job logs via `/api/jobs/{id}` endpoint. Common issues:
- GECKO-A executable not found
- Insufficient disk space
- Invalid VOC name

### Viewing Downloaded 3D Files (MOL/SDF/XYZ)
Use molecular visualization software:
- **Avogadro** (free, recommended): https://avogadro.cc/
- **PyMOL**: https://pymol.org/
- **UCSF ChimeraX**: https://www.rbvi.ucsf.edu/chimera/
- **Online**: MolView.org (paste file content)

## License

This project is provided for research and educational purposes.
