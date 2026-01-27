# GECKO-A Web Interface - Claude Code Instructions

## Project Overview

GECKO-A is an atmospheric chemistry web application that generates chemical mechanisms and runs box model simulations for Volatile Organic Compounds (VOCs).

**Tech Stack:**
- **Backend:** FastAPI (Python 3.9+)
- **Frontend:** Vanilla JS with Cytoscape.js, 3Dmol.js
- **Scientific:** GECKO-A & BOXMODEL4GECKO (Fortran 90)
- **Visualization:** RDKit, Graphviz, Matplotlib

---

## Key Files to Read First

### Understanding the Application
```
gecko_web/main.py              # FastAPI backend - ALL endpoints (3000+ lines)
gecko_web/templates/index.html # Single-page application UI
gecko_web/static/js/app.js     # Frontend logic (48KB)
```

### Configuration & Deployment
```
fly.toml                       # Fly.io deployment config
docker/Dockerfile              # Container build instructions
docker/docker-compose.yml      # Local Docker setup
FLY_DEPLOYMENT.md              # Step-by-step cloud deployment guide
DEPLOYMENT_INSTRUCTIONS.md     # Full deployment documentation
```

### Scientific Core (Fortran)
```
docker/gecko_source/           # GECKO-A mechanism generator
docker/gecko_source/LIB/       # Fortran 90 source modules
docker/gecko_source/INPUT/     # Configuration templates (gecko.nml)
docker/boxmodel_source/        # Box model solver
```

### Chemical Database
```
gecko_web/chemdata/compound_database.py  # 150+ VOC compounds
gecko_web/chemdata/voc_categories.py     # VOC categorization
gecko_web/chemdata/reaction_data.py      # Kinetics data
```

---

## Common Tasks

### 1. Add a New VOC Compound
Edit: `gecko_web/chemdata/compound_database.py`
```python
"compound_name": {
    "smiles": "SMILES_STRING",
    "gecko_formula": "GECKO_NOTATION",
    "molecular_weight": 123.45,
    "category": "alkane|alkene|aromatic|terpene|oxygenated"
}
```

### 2. Modify API Endpoints
Edit: `gecko_web/main.py`
- Job endpoints: Lines 400-600
- Compound endpoints: Lines 600-800
- Workflow endpoints: Lines 800-1000

### 3. Update Frontend UI
Edit: `gecko_web/templates/index.html` (HTML structure)
Edit: `gecko_web/static/js/app.js` (JavaScript logic)
Edit: `gecko_web/static/css/styles.css` (Styling)

### 4. Modify Docker Build
Edit: `docker/Dockerfile`
- Python dependencies: Line 18-19
- Fortran compilation: Lines 26-37

### 5. Change Fly.io Settings
Edit: `fly.toml`
- Memory: `memory_mb` under `[[vm]]`
- Region: `primary_region`

---

## Running Locally

### Option 1: macOS (UI Only)
```bash
./start_gecko.command
# Opens http://localhost:8000
```

### Option 2: Docker (Full Functionality)
```bash
docker-compose -f docker/docker-compose.yml up
# Opens http://localhost:8000
```

### Option 3: Python Dev Mode
```bash
source venv/bin/activate
pip install -r gecko_web/requirements.txt
uvicorn gecko_web.main:app --reload --port 8000
```

---

## Deploying to Fly.io

### First Time
```bash
fly auth login
fly launch --no-deploy
fly volumes create gecko_data --region iad --size 10
fly deploy
```

### Subsequent Deploys
```bash
git add . && git commit -m "message" && git push
# GitHub Actions auto-deploys, OR:
fly deploy
```

---

## Architecture

```
User Browser
     │
     ▼
┌─────────────────────────────────────┐
│  FastAPI (gecko_web/main.py)        │
│  - REST API endpoints               │
│  - Job management                   │
│  - Static file serving              │
└─────────────────────────────────────┘
     │
     ├──► RDKit (molecular structures)
     ├──► Graphviz (pathway diagrams)
     ├──► Matplotlib (plots)
     │
     ▼
┌─────────────────────────────────────┐
│  GECKO-A (Fortran)                  │
│  - Mechanism generation             │
│  - Reaction pathway expansion       │
└─────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────┐
│  BOXMODEL4GECKO (Fortran)           │
│  - Kinetic simulations              │
│  - Time-series output               │
└─────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────┐
│  /data/output/{job_id}/             │
│  - .mech, .kpp files                │
│  - .nc (NetCDF) results             │
│  - .png diagrams                    │
└─────────────────────────────────────┘
```

---

## API Endpoints Quick Reference

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/` | GET | Web UI |
| `/api/jobs` | GET | List all jobs |
| `/api/jobs` | POST | Create new job |
| `/api/jobs/{id}` | GET | Job status |
| `/api/jobs/{id}/results` | GET | Job results |
| `/api/compounds` | GET | All compounds |
| `/api/compounds/search/{q}` | GET | Search |
| `/api/kinetics/{name}` | GET | Rate constants |
| `/api/environment` | GET | System status |

---

## Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `DATA_DIR` | `/data` | Job output storage |
| `GECKO_SOURCE_DIR` | `/app/gecko_source` | GECKO-A location |
| `BOXMODEL_SOURCE_DIR` | `/app/boxmodel_source` | Box model location |
| `PORT` | `8000` | Server port |

---

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test
pytest tests/test_api.py -v

# With coverage
pytest tests/ --cov=gecko_web
```

---

## Troubleshooting

### Build fails on Fly.io
- Check `fly logs` for errors
- Ensure Dockerfile has all system dependencies
- Fortran compilation takes 45-90 minutes (normal)

### Simulations fail
- Check memory: Complex VOCs need 4GB+ RAM
- Check disk: Need 10GB+ for job storage
- Check `/api/environment` for system status

### Frontend not loading
- Check browser console for JS errors
- Verify CDN libraries loading (Cytoscape, 3Dmol)
- Check `/api/jobs` returns data

---

## Files NOT to Modify

- `docker/gecko_source/LIB/*.f90` - Core Fortran (unless you know Fortran)
- `docker/boxmodel_source/LIBSRC/*` - Box model core
- Compiled files (`*.o`, `*.mod`) - Auto-generated

---

## Version History

See `CHANGELOG.md` for full history.

Current: **v3.0.6** (2026-01-27)
