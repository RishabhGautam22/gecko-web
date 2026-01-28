# What Gets Uploaded to GitHub

## Files INCLUDED (Will Be on GitHub & Fly.io)

### Python Application (~1.5 MB)
```
gecko_web/
├── main.py                    # FastAPI backend
├── mechanism_diagram.py       # Diagram generation
├── pathway_visualizer.py      # RDKit visualization
├── enhanced_visualizer.py     # Advanced visualization
├── reaction_tree.py           # Reaction analysis
├── postprocessing.py          # Data processing
├── mass_balance.py            # Atom conservation
├── requirements.txt           # Python dependencies
├── chemdata/
│   ├── compound_database.py   # 150+ VOC compounds
│   ├── voc_categories.py      # VOC categories
│   └── reaction_data.py       # Kinetics data
├── templates/
│   └── index.html             # Web UI
└── static/
    ├── css/styles.css
    └── js/*.js                # Frontend scripts
```

### Fortran Source Code (~147 MB)
```
docker/gecko_source/           # GECKO-A (11 MB)
├── LIB/*.f90                  # 80+ Fortran modules ✅ INCLUDED
├── DATA/*.dat                 # Reaction databases ✅ INCLUDED
├── INPUT/                     # Config templates ✅ INCLUDED
│   ├── gecko.nml
│   └── cheminput.dat
├── RUN/                       # Run scripts ✅ INCLUDED
└── OBJ/makefile              # Build config ✅ INCLUDED

docker/boxmodel_source/        # Box Model (136 MB)
├── LIBSRC/*.f90              # Source modules ✅ INCLUDED
├── PROG/*.f90                # Main programs ✅ INCLUDED
├── SCRIPTS/                  # Build scripts ✅ INCLUDED
└── build.sh                  # Compilation script ✅ INCLUDED
```

### Configuration & Deployment
```
fly.toml                       # Fly.io config ✅
docker/Dockerfile              # Container build ✅
docker/docker-compose.yml      # Local Docker ✅
docker/docker_arch.gnu         # NetCDF config ✅
.github/workflows/*.yml        # CI/CD pipelines ✅
```

### Documentation
```
README.md                      ✅
DOCUMENTATION.md               ✅
CHANGELOG.md                   ✅
DEPLOYMENT_INSTRUCTIONS.md     ✅
FLY_DEPLOYMENT.md              ✅
CLAUDE.md                      ✅ (VS Code instructions)
```

### Scripts
```
scripts/start.sh               # Unix launcher ✅
scripts/start.bat              # Windows launcher ✅
scripts/start.ps1              # PowerShell launcher ✅
start_gecko.command            # macOS launcher ✅
```

### Sample Data
```
samples/                       # Example input files ✅
```

---

## Files EXCLUDED (Not Uploaded)

### Large/Generated Files
```
.venv/                         # Virtual environment (465 MB) ❌
venv/                          # Alternative venv ❌
Gecko Windows Installation/    # Windows installer (645 MB) ❌
TestGeko/                      # Test data (59 MB) ❌
test_environment/              # Test infrastructure ❌
test_vis_bench/                # Visualization benchmarks ❌
```

### Build Artifacts (Rebuilt in Docker)
```
docker/gecko_source/OBJ/*.o    # Compiled objects ❌
docker/gecko_source/OBJ/*.mod  # Fortran modules ❌
docker/gecko_source/OBJ/cm     # GECKO executable ❌
docker/boxmodel_source/OBJ/*   # Compiled objects ❌
docker/boxmodel_source/PROG/boxmod  # Box model executable ❌
```

### User Data
```
data/output/*                  # Job results ❌
data/mechanisms/*              # Cached mechanisms ❌
data/archives/*                # Archived jobs ❌
```

### System/IDE Files
```
.DS_Store                      ❌
__pycache__/                   ❌
.pytest_cache/                 ❌
.env                           # Secrets ❌
```

---

## Estimated Repository Size

| Component | Size |
|-----------|------|
| gecko_web/ | ~1.5 MB |
| docker/gecko_source/ | ~11 MB |
| docker/boxmodel_source/ | ~30 MB (source only, no CHEMDAT) |
| Documentation | ~100 KB |
| Configuration | ~50 KB |
| **Total** | **~45-50 MB** |

Note: The boxmodel_source CHEMDAT directory (226 mechanism directories) is excluded by its own .gitignore to keep repo size manageable.

---

## Verification Commands

```bash
# Check what will be committed
git status

# Check total size of tracked files
git ls-files | xargs du -ch | tail -1

# Verify Fortran sources are tracked
git ls-files | grep "\.f90$" | wc -l
# Should show ~492 files

# Verify key files exist
git ls-files docker/gecko_source/LIB/
git ls-files docker/gecko_source/DATA/
```

---

## What Happens on Fly.io

When you deploy to Fly.io:

1. **GitHub → Fly.io**: All tracked files are sent
2. **Docker Build**:
   - Installs Python dependencies from `requirements.txt`
   - Compiles Fortran from source (`docker/gecko_source/LIB/*.f90`)
   - Compiles Box Model from source (`docker/boxmodel_source/`)
3. **Runtime**: Full GECKO-A functionality available

The Fortran executables are **built fresh** during Docker build, not uploaded pre-compiled. This ensures compatibility across platforms.
