# Changelog

All notable changes to the GECKO-A Web Interface are documented in this file.

**Author: Deeksha Sharma**

## [3.0.10] - 2026-01-29

### Chemical Structure Visualization Improvements

#### Fixed Terpene Oxidation Product Ring Structures (CRITICAL)
- **File**: `gecko_web/reaction_tree.py`
- Fixed **pinic acid**, **pinonic acid**, **pinonaldehyde** to show correct **cyclobutane (4-membered) rings**
- Previously displayed incorrect cyclohexane (6-membered) rings due to fallback SMILES generation
- Updated `KNOWN_SPECIES` dictionary with PubChem-verified SMILES:
  - Pinic acid: C9H14O4 (CID 10131) - cyclobutane dicarboxylic acid
  - Pinonic acid: C10H16O3 (CID 10130) - cyclobutane ketoacid
  - Pinonaldehyde: C10H16O2 (CID 17616) - cyclobutane aldehyde
- Added GECKO species codes (TA9000, TA8000, TA7000) for proper lookup
- Improved `_generate_fallback_smiles()` to detect ring sizes from GECKO notation

#### Human-Readable Compound Names in Diagrams
- **File**: `gecko_web/reaction_tree.py`
- Added comprehensive `GECKO_CODE_TO_NAME` mapping dictionary (150+ entries)
- New `get_compound_name()` function converts GECKO codes to human-readable names:
  - `ISOPRE` → "Isoprene"
  - `2U5000` → "Isoprene-OH adduct"
  - `MVK` → "Methyl vinyl ketone"
  - `TA9000` → "Pinic acid"
  - `R07000` → "Toluene"

- **File**: `gecko_web/mechanism_diagram.py`
- Updated `generate_static_diagram()` to use `get_compound_name()`
- Title and labels now show human-readable names

- **File**: `gecko_web/pathway_visualizer.py`
- Updated `from_reaction_tree()` to convert GECKO codes to readable names
- Imported `get_compound_name` from reaction_tree

#### Added Compound Labels to Diagrams
- **File**: `gecko_web/pathway_visualizer.py`
- Added compound name (dark blue) and molecular formula (gray) below each structure
- New `_add_labels_to_mol_image()` method using PIL for text overlay
- Labels are properly centered and sized for readability

- **File**: `gecko_web/mechanism_diagram.py`
- Updated `_smiles_to_mol_image()` to accept name and formula parameters
- Labels displayed consistently with pathway diagrams

#### Comprehensive 130-Compound Validation
- **File**: `tests/comprehensive/generate_all_130_compounds.py` (NEW)
- Generates labeled PNG images for ALL 130 dropdown compounds
- Validates SMILES, molecular formula, ring structures, and molecular weight
- Outputs JSON validation report with detailed compound data

#### Added Missing Compounds
- **File**: `gecko_web/chemdata/compound_database.py`
- Added **nerolidol** (C15H26O) - sesquiterpene alcohol
- Added **butyl_acetate** (C6H12O2) - ester
- Added **neopentane** (C5H12) - branched alkane

### Validation Results
- All 130 compounds generate valid PNG images
- Critical terpene oxidation products confirmed with cyclobutane rings:
  - pinic_acid: rings=[4], C9H14O4 ✓
  - pinonic_acid: rings=[4], C10H16O3 ✓
  - pinonaldehyde: rings=[4], C10H16O2 ✓
  - nopinone: rings=[4, 6, 6], C9H14O ✓

### Known Limitations
- **Box model failures**: Some compounds (e.g., nerolidol, butyl_acetate) may fail box model simulations due to GECKO-A's internal species dictionary limitations. The mechanism generator still produces valid outputs for these compounds.

---

## [3.0.9] - 2026-01-28

### Scientific Audit Fixes

PhD-level scientific review of physics/chemistry core with 6 critical fixes:

#### SOA Yield Placeholder Warnings (CRITICAL)
- **File**: `gecko_web/main.py` (lines 2193-2236)
- SOA yield comparison function now includes prominent warnings
- Added `is_placeholder: true` flag and `data_type: "literature_estimate"`
- Literature sources cited for each VOC category (Ng et al., Griffin et al., etc.)
- Users now see explicit warnings that yields are NOT from simulation

#### Vapor Pressure Fallback Warnings (CRITICAL)
- **File**: `gecko_web/postprocessing.py`
- Added critical warning when `pvap.nan.dat` file is missing
- Fallback to MW-based estimation now explicitly flagged
- `data_quality` object now includes `pvap_source` and `pvap_accuracy_warning`

#### Synthetic Data Plot Watermarks
- **File**: `gecko_web/postprocessing.py`
- New `_add_synthetic_watermark()` function for all plot generation
- Watermarks applied when data is synthetic, surrogate, or has critical Pvap quality
- Red "SYNTHETIC DATA", orange "SURROGATE DATA", yellow "ESTIMATED PVAP"

#### Frontend Data Quality Warnings
- **File**: `gecko_web/static/js/app.js`
- Enhanced `renderDataQuality()` with severity levels (info/warning/critical)
- Vapor pressure quality indicator with color-coded warnings
- SOA yield placeholder warnings in comparison results
- **File**: `gecko_web/static/css/styles.css`
- New CSS classes: `.severity-critical`, `.severity-warning`, `.soa-yield-warning`

#### Configurable Mass Balance Tolerances
- **File**: `gecko_web/mass_balance.py`
- `MassBalanceChecker` now accepts configurable tolerances via constructor
- Added factory methods: `with_strict_tolerances()`, `with_relaxed_tolerances()`
- Default, strict (0.001), and relaxed (0.1) tolerance presets

#### Arrhenius Validation Tests
- **File**: `tests/test_arrhenius_nist.py` (NEW)
- 20 comprehensive tests validating against NIST/JPL reference data
- Tests for OH + CH4, OH + C2H6, OH + C3H8, OH + CO, O3 + alkenes, etc.
- Temperature dependence validation and uncertainty bounds checking
- Troe pressure-dependent reaction tests

### Technical Details
- All 171 tests pass (154 passed, 17 skipped)
- No breaking changes to existing functionality
- Follows "Surgical Changes" and "Simplicity First" principles

---

## [3.0.8] - 2026-01-28

### Baseline Manager & Data Quality Appendix
- New Baseline Manager for automatic surrogate data fallback
- Data Quality Appendix PDF generated for every simulation
- Enhanced data provenance tracking

---

## [3.0.7] - 2026-01-28

### Code Cleanup & Optimization

#### Removed Obsolete Files (14 files deleted)
- **Debug/Test Scripts**: Removed `regenerate_diagrams.py`, `regenerate_trees_v2.py`, `test_parser.py`, `test_pipeline_2.py`, `debug_parser_2.py`, `inspect_json.py`, `inspect_keys.py`, `fix_current_job.py`, `test.f90`
- **Unused Deployment Configs**: Removed `railway.toml`, `render.yaml` (project standardized on Fly.io)
- **Stale Artifacts**: Removed `server.log`, `test_bin`, `test_smiles.html`
- **Unused JavaScript Module**: Removed `gecko_web/static/js/structure3d.js` (never imported)

#### Cleaned Python Imports (11 files updated)
- Removed unused imports across all Python modules including:
  - `numpy`, `pandas` from `combined_workflow.py`
  - `json`, `Path`, `Descriptors` from `compound_database.py`
  - `io`, `math`, `networkx`, `PIL` from `enhanced_visualizer.py`
  - `glob`, `tempfile`, `contextmanager` from `main.py`
  - Many more unused typing imports and RDKit submodules

#### Cleaned JavaScript API (api.js)
- Removed 11 unused API functions:
  - `getCompounds()`, `getCompoundCategories()`, `getCompound()`, `searchCompounds()`, `getCompoundsByCategory()`
  - `getKinetics()`, `getAtmosphericLifetime()`, `createCombinedWorkflow()`
  - `getJob3DStructures()`, `pollJobStatus()`, `checkEnvironment()`

#### Data Cleanup
- Reset `data/job_state.json` (had stale paths to old directory)
- Cleared old job outputs from `data/output/`

### Technical Details
- All 142 tests pass (125 passed, 17 skipped)
- All Python modules compile without errors
- All JavaScript files pass syntax validation
- No breaking changes to existing functionality

---

## [3.0.6] - 2026-01-27

### Enhancements
- **Consistent Pathway Diagrams**: 
    - Implemented `enhanced_visualizer` for both Generator and Box Model jobs.
    - Diagrams now feature full-color CPK rendering and functional group labels across all workflow steps.
- **Improved 3D Visualization**:
    - **CPK Standards**: Carbon atoms are now rendered in dark charcoal (#1A1A1A) to match standard chemical visualization conventions.
    - **Descriptive Naming**: The 3D viewer now displays the chemical code alongside detected functional groups (e.g., "U5000 [Nitrate, Alcohol]") instead of opaque codes.
- **Robust Download Links**: Fixed an issue where the "Download PNG" button served legacy black-and-white diagrams or broken links.

### Bug Fixes
- **Double Extension Fix**: Resolved a bug where diagrams were saved as `.png.png`, preventing them from appearing in the UI.
- **Functional Group Preservation**: Fixed a regression where functional group labels were being overwritten by empty lists in the visualizer.

## [3.0.5] - 2026-01-27

### Critical Fixes
- **Visualizer Structure Correction**: Fixed issue where branched alkanes (Isopentane) were reverting to linear structures in generated diagrams. The visualizer now strictly respects the SMILES generated by the `reaction_tree` parser.
- **Graphviz Syntax Error**: Resolved a crash in diagram generation caused by duplicate attribute definitions (`splines=...`) passed to the Graphviz engine.
- **SMILES Priority Logic**: Updated `pathway_visualizer.py` to prioritize explicit input SMILES over static dictionary lookups for ambiguous GECKO codes like `C05000`.

## [3.0.4] - 2026-01-24

### Breaking Changes

#### Removed Combined Workflow
- **Combined Workflow has been removed**: The complex combined workflow feature was causing failures and has been completely removed
- **Use Generator then Box Model separately**: Run Generator first to create mechanism, then Box Model to simulate
- **Simplified UI**: Job type selector now shows only Generator, Box Model, and Comparison modes

### Bug Fixes

#### Fixed Isopentane/Branched Alkane Structure Display
- **Added GECKO formula to SMILES mappings**: Added direct mappings for GECKO formulas like `CH3CH(CH3)CH2CH3`
- **Extended special_mappings**: Added isopentane, neopentane, isobutane, cyclohexane, cyclopentane
- **Improved root node detection**: Better matching for branched compounds

#### Fixed Pathway Diagram Chaos
- **Reduced max_nodes default**: Changed from 200 to 50 nodes for cleaner diagrams
- **Improved graphviz settings**: Changed splines from 'polyline' to 'ortho' for cleaner edges
- **Added edge concentration**: Merge edges going to same node with `concentrate='true'`
- **Increased spacing**: Better nodesep (0.8) and ranksep (1.5) for readability

### Technical Details

- Modified: `gecko_web/main.py` - Removed combined workflow import, endpoint, and request model
- Modified: `gecko_web/templates/index.html` - Removed Combined Workflow option from dropdown
- Modified: `gecko_web/static/js/app.js` - Removed combined workflow handling code
- Modified: `gecko_web/reaction_tree.py` - Added GECKO formula mappings, reduced max_nodes
- Modified: `gecko_web/pathway_visualizer.py` - Improved graphviz layout settings

---

## [3.0.3] - 2026-01-24

### Critical Bug Fixes

#### Simplified Combined Workflow (FIXED)
- **Completely rewrote combined workflow**: Removed complex multi-stage workflow orchestration
- **Simple sequential execution**: Generator runs first, then Box Model when it completes
- **Single output directory**: All files go to one location for easy access
- **Better error handling**: Clear error messages at each stage

#### 3D Structure Generation for Radical Species (FIXED)
- **Improved SMILES sanitization**: Better handling of radical oxygen `[O]`, peroxy radicals `O[O]`, and nitrogen radicals
- **Extended pattern matching**: Handles terminal radicals, internal radicals, and carbon radicals
- **More robust embedding**: Multiple fallback strategies for difficult structures

#### Fixed GECKO Formulas for Cyclic Compounds (FIXED)
- **Fixed cyclohexane**: Changed from `C1H2CH2CH2CH2CH2CH21` to `C1H2CH2CH2CH2CH2C1H2`
- **Fixed cyclopentane**: Changed from `C1H2CH2CH2CH2CH21` to `C1H2CH2CH2CH2C1H2`
- **Fixed cyclohexene**: Changed from `C1H2CH2CH2CdH=CdHCH21` to `C1H2CH2CH2CdH=CdHC1H2`
- **Fixed tetrahydrofuran**: Changed from `C1H2CH2CH2OCH21` to `C1H2CH2CH2OC1H2`
- **Fixed 3-carene**: Now uses correct GECKO bicyclic notation from fixedname.dat
- **Fixed many terpenes**: alpha_terpinene, gamma_terpinene, terpinolene, alpha_terpineol, etc.
- **Fixed sesquiterpenes**: beta_caryophyllene, alpha_humulene
- **Fixed oxidation products**: pinonaldehyde, pinonic_acid, pinic_acid

### Technical Details

- Modified: `gecko_web/main.py` - Simplified `run_combined_workflow()`, improved `_sanitize_smiles_for_3d()`
- Modified: `gecko_web/chemdata/compound_database.py` - Fixed 15+ cyclic compound GECKO formulas

---

## [3.0.2] - 2026-01-23

### Critical Bug Fixes

#### Combined Workflow Failures (FIXED)
- **Fixed mechanism directory passing in combined workflow**: The `_stage_box_model` method was passing incorrect directory to `run_box_model`, causing box model preparation to fail because it couldn't find mechanism files
- **Added new `run_box_model_with_mechanism_dir()` function**: Accepts separate `mechanism_dir` and `output_dir` parameters for proper combined workflow execution
- **Improved error handling**: Better logging and error messages for debugging mechanism import failures

#### Box Model Simulation Failures (FIXED)
- **Fixed `get_gecko_input()` function**: Now properly checks the compound database for GECKO formulas instead of only using the hardcoded VOC_MAPPING
- **Extended compound support**: Compounds like isobutane, cyclohexane, and many others now work correctly
- **Added alternate name format support**: Handles hyphens, underscores, and spaces in compound names

#### 3D Model Positioning (FIXED)
- **Fixed 3D viewer container positioning**: Added CSS styles to properly contain the $3Dmol canvas within its dedicated container
- **Prevented canvas escape**: Set `position: relative` on container and `position: absolute` on canvas to keep 3D viewer in place
- **Improved JavaScript initialization**: Better container dimension handling before viewer creation

#### 3D Model Loading Failures (FIXED)
- **Added SMILES sanitization for 3D embedding**: Converts radical notation to embeddable forms for visualization
- **Implemented multiple embedding strategies**:
  1. Standard embedding
  2. Random coordinates
  3. ETKDG with increased iterations
  4. Sanitized SMILES with random coordinates
  5. 2D coordinates projected to 3D (last resort)
- **Improved error messages**: Clear feedback about why certain structures (radicals, unusual valence) cannot be embedded

### Technical Details

- Modified: `gecko_web/combined_workflow.py` - Updated `_stage_box_model` to use new function
- Modified: `gecko_web/main.py` - Added `run_box_model_with_mechanism_dir()` and `_sanitize_smiles_for_3d()` functions, improved `get_gecko_input()`
- Modified: `gecko_web/static/css/styles.css` - Added 3D viewer container styles
- Modified: `gecko_web/static/js/app.js` - Improved `load3DStructure()` function

---

## [3.0.1] - 2026-01-18

### Bug Fixes

#### Combined Workflow Fixes
- **Fixed WorkflowConfig attributes**: Added missing `output_base`, `run_generator`, `run_boxmodel`, `run_postprocessing`, `verify_mass_balance`, `generator_options`, and `boxmodel_options` fields
- **Fixed CombinedWorkflow initialization**: Now accepts `output_dir` from config or parameter
- **Fixed OutputStructure**: Uses base directory directly instead of creating nested workflow_id folder
- **Fixed WorkflowResult**: Added `success` and `error` properties for proper status checking
- **Fixed JSON serialization**: Path objects now properly converted for JSON output

#### Enhanced Visualizer Fixes
- **Fixed graphviz type hints**: Changed `graphviz.Digraph` to string quotes to prevent NameError when graphviz not installed
- **Added `generate_diagram()` method**: API-compatible wrapper for the `generate()` method
- **Added DiagramStyle attributes**: `show_branching_ratios` and `color_by_reaction_type` fields

#### 3D Structure Viewer Fixes
- **Fixed API response handling**: JavaScript now correctly uses `response.data` instead of `response.structure`
- **Added SMILES pass-through**: 3D endpoint now accepts SMILES strings directly
- **Improved error messages**: Clear feedback when RDKit is not available

#### UI Improvements
- **Added VOC Comparison Mode**: Multi-select dropdown for comparing 2-10 VOCs
- **Added Mechanism Reduction UI**: DRGEP, PFA, Lumping, Sensitivity-based methods with configurable parameters
- **Added 3D Structure Viewer**: Interactive 3Dmol.js viewer with species dropdown
- **Added Layout Algorithm selector**: Hierarchical, Radial, Force-Directed options for pathway diagrams
- **Added PDF Report downloads**: Mass Balance Report, Mechanism Summary (JSON)
- **Added Mass Balance display**: Shows C/H/O/N balance percentages with status

### Dependencies
- **rdkit**: Now properly required for 3D structure generation
- **graphviz**: Both system package and Python bindings required for diagram regeneration

---

## [3.0.0] - 2026-01-18

### Major Scientific Overhaul

This release implements 15 comprehensive improvements to scientific accuracy, visualization, and usability based on an independent code audit.

### New Modules

#### chemdata/ Package
- **compound_database.py**: 150+ VOC compounds with verified SMILES, molecular properties, and reaction kinetics
  - Validated against NIST, PubChem, and EPA CompTox
  - Includes CAS numbers, InChI identifiers, and literature references
  - Categories: alkanes, alkenes, aromatics, terpenes, oxygenated VOCs

- **voc_categories.py**: Hierarchical VOC categorization
  - Alkanes (linear, branched, cyclic)
  - Alkenes (simple, dienes, cyclic)
  - Aromatics (BTEX, PAHs)
  - Terpenes (monoterpenes, sesquiterpenes)
  - Oxygenated (alcohols, aldehydes, ketones, acids)

- **reaction_data.py**: Atmospheric reaction kinetics database
  - Arrhenius parameters for OH, O3, NO3 reactions
  - Troe pressure-dependent parameters
  - Temperature-dependent rate calculations
  - IUPAC/JPL recommended values

#### Core Modules
- **mass_balance.py**: Mass balance verification for chemical mechanisms
  - Atom counting (C, H, O, N, S)
  - Reaction balance checking
  - Carbon tracking through pathways
  - JSON report generation

- **combined_workflow.py**: Unified Generator + Box Model workflow
  - Stage-based execution with progress tracking
  - Structured output directories (mechanism/, simulation/, analysis/, plots/, reports/)
  - Error handling and recovery
  - PDF report generation with ReportLab

- **enhanced_visualizer.py**: Publication-quality pathway visualization
  - RDKit-based molecular rendering
  - Reaction type color coding:
    - OH addition: blue (#3498db)
    - OH abstraction: green (#2ecc71)
    - O3 ozonolysis: orange (#e67e22)
    - NO3 reactions: purple (#9b59b6)
    - Photolysis: yellow (#f1c40f)
    - Isomerization: cyan (#1abc9c)
  - Multiple layout algorithms (hierarchical, radial, force-directed)
  - Branching ratios on edges
  - Variable legend (X = -OOH, -ONO2, -OH)
  - Stereochemistry detection (E/Z isomers)

#### Frontend Modules
- **static/js/api.js**: Extended API client with all new endpoints
- **static/js/structure3d.js**: 3D molecular visualization with 3Dmol.js

### New API Endpoints

#### Compound Database
- `GET /api/compounds` - List all compounds with properties
- `GET /api/compounds/categories` - Get VOC category hierarchy
- `GET /api/compounds/{name}` - Get detailed compound info
- `GET /api/compounds/search/{query}` - Search compounds
- `GET /api/compounds/category/{cat}` - Get compounds by category

#### Reaction Kinetics
- `GET /api/kinetics/{name}` - Get rate constants at temperature
- `GET /api/kinetics/{name}/lifetime` - Calculate atmospheric lifetime

#### Advanced Workflows
- `POST /api/workflow/combined` - Combined Generator + Box Model
- `POST /api/workflow/comparison` - Compare multiple VOCs
- `POST /api/mechanism/reduce` - Reduce mechanism size
- `GET /api/jobs/{id}/mass-balance` - Verify mass balance
- `POST /api/jobs/{id}/visualize` - Generate enhanced diagram

#### 3D Structures
- `GET /api/structure/{name}/3d` - Get 3D structure (MOL/SDF/XYZ/PDB)
- `GET /api/jobs/{id}/structures/3d` - Get all structures for job

### New Pydantic Models

- `ExtendedGeckoOptions` - Full GECKO-A generator options
  - Vapor pressure method selection (Nannoolal, SIMPOL, etc.)
  - Reaction channel toggles (OH, O3, NO3, photolysis)
  - SAR method selection
  - Criegee chemistry options
  - Autoxidation parameters

- `ExtendedBoxModelOptions` - Full box model configuration
  - Environmental conditions (T, P, RH, location, time)
  - Initial concentrations (VOC, O3, NOx, etc.)
  - Aerosol parameters (seed, density, activity coefficient)
  - Dilution and deposition settings

- `CombinedWorkflowRequest` - Unified workflow configuration
- `VOCComparisonRequest` - Multi-VOC comparison settings
- `MechanismReductionRequest` - Reduction parameters (DRGEP, lumping, PFA)

### Mechanism Reduction Methods

Implemented four reduction algorithms:
1. **DRGEP** (Directed Relation Graph with Error Propagation)
   - Based on Pepiot-Desjardins & Pitsch (2008)
   - Path-based importance analysis
   - Configurable error threshold

2. **Lumping** - Species grouping by structure
   - Groups by carbon number and functional groups
   - Preserves radical/SOA-precursor species

3. **PFA** (Path Flux Analysis)
   - Based on Sun et al. (2010)
   - Flux-weighted species importance

4. **Sensitivity-based** - Connectivity analysis
   - Simple reaction graph connectivity
   - Configurable minimum species

### Improvements

- **requirements.txt**: Updated with version constraints and new dependencies
- **start_gecko.command**: Updated for v3.0.0 with new feature descriptions
- Comprehensive error handling throughout new modules
- Type hints and docstrings for all new code

---

## [2.0.6] - 2026-01-18

### Bug Fixes

#### Pathway Diagram Rendering (CRITICAL)
- **Fixed duplicate edges in reaction tree**: Removed duplicate edges that caused graphviz layout distortion
  - Previously: 234 edges with many duplicates (e.g., 33 edges between same node pairs)
  - Cause: `_build_tree_bfs` in reaction_tree.py added edges for every reaction without deduplication
  - Result: Molecule images were compressed/cropped, making ethylbenzene appear as benzene
  - Now: 27 unique edges, molecules render correctly with full substituent groups visible

- **Root cause**: When graphviz renders many overlapping arrows between the same nodes, the layout algorithm distorts node sizes and positions, causing molecule PNG images to be scaled incorrectly

- **Fix applied**: Added `seen_edges` set to track unique (from, to) pairs and aggregate yields from duplicate reactions

### Verification
- Ethylbenzene root molecule now correctly shows full ethyl group (-CH2CH3)
- All aromatic substituents render with correct carbon chain lengths
- Pathway diagrams maintain proper molecular structure proportions

### Windows Native Compilation Support
- **Complete GECKO-A source code**: All 81 Fortran 90 source files included
- **Complete Box Model source code**: All source files for BOXMODEL4GECKO
- **BUILD_GECKO.bat**: Automated build script for Windows with MinGW-w64
- **Makefile.windows**: Windows-compatible makefile (removes macOS -no_pie flag)
- **gecko.bat**: Windows batch wrapper script for GECKO-A execution
- **Pre-computed CHEMDAT**: Mechanism files for isoprene, alpha-pinene, toluene, ethane
- **Improved path handling**: Proper Windows path detection and temp directory usage
- **Better error messages**: Clear instructions when executables are missing

### Build Requirements (for native compilation)
- MSYS2 with MinGW-w64 (gfortran compiler)
- NetCDF libraries (for Box Model)
- See README.txt in Windows Installation folder for detailed instructions

---

## [2.0.5] - 2026-01-18

### Major Improvements

#### SMILES Generation Overhaul
- **Expanded KNOWN_SPECIES database**: From ~100 to 500+ validated entries
  - Complete coverage of C1-C15 linear alkanes and oxidation products
  - Terpenes (monoterpenes, sesquiterpenes) and their derivatives
  - Aromatic compounds including all xylene isomers and toluene products
  - Common atmospheric radicals and intermediates

- **Unified SMILES conversion**: Single authoritative `gecko_to_smiles()` function in reaction_tree.py
  - All modules now delegate to this function
  - Eliminates duplicate/conflicting conversion code
  - Consistent SMILES output across the entire application

- **Fixed aromatic pattern parser**: Proper regex for GECKO's benzene ring notation
  - Handles `c1HcHcHcHcHc1` patterns correctly
  - Supports substituted aromatics

- **RDKit validation layer**: Optional SMILES validation and canonicalization
  - Validates all generated SMILES when RDKit is available
  - Falls back gracefully when RDKit is not installed

#### Frontend Modernization
- **Refactored index.html**: Reduced from 1127 lines to 290 lines
  - Extracted JavaScript into modular files:
    - `api.js`: API communication
    - `reaction-tree.js`: Cytoscape visualization
    - `app.js`: Main application logic
  - Extracted CSS into `styles.css`

- **Enhanced UI design**: Modern CSS with CSS variables
  - Consistent color scheme and spacing
  - Improved responsive design
  - Loading animations
  - Better status badge styling

- **Added author attribution**: Header and footer with version info

#### Documentation
- **Combined documentation**: All .md files merged into DOCUMENTATION.md
- **Removed false claims**: Updated README to reflect actual capabilities
- **Added comprehensive SMILES unit tests**: tests/test_smiles.py

### Bug Fixes
- Fixed dead code removal in mechanism_diagram.py
- Fixed import statements for unified SMILES function

---

## [2.0.4] - 2026-01-16

### New Features

#### Publication-Quality Pathway Diagrams
- **New `pathway_visualizer.py` module**: RDKit + Graphviz-based diagram generator
  - Produces GECKO-A publication-style oxidation mechanism diagrams
  - Hierarchical tree layout with molecular structure images at each node
  - Branching ratios displayed on edges
  - Support for PNG and SVG output formats
  - Automatic SMILES cleaning for RDKit compatibility (radical dots, nitrate groups, cumulated double bonds)

### Bug Fixes

#### Cytoscape Reaction Tree Visualization (CRITICAL)
- **Fixed edge validation**: Edges now only reference nodes that exist in the graph
- **Fixed edge ordering**: Nodes are now added before edges reference them

#### Aromatic SMILES Generation (CRITICAL)
- **Fixed aromatic ring SMILES conversion**: Toluene and derivatives now render correctly
- Added 15+ toluene oxidation products to KNOWN_SPECIES

### Improvements
- Comprehensive frontend logging for Cytoscape initialization
- Two-stage layout execution for improved reliability

---

## [2.0.3] - 2026-01-16

### Improvements
- **Increased DPI to 300**: All plots now publication-ready
- **Scientific Validation Tests**: 19 comprehensive tests added

### Bug Fixes
- **Fixed VOC Root Detection**: Now reads from `listprimary.dat`

---

## [2.0.2] - 2025-01-16

### Critical Data Fixes
- **netCDF Reading**: Fixed concentration data reading
- **Vapor Pressure**: Fixed interpretation (was off by ~10 orders of magnitude)
- **C* Calculation**: Now matches Pankow (1994)
- **Dictionary Parsing**: Fixed carbon counts and molecular weights

---

## [2.0.1] - 2025-01-15

### Bug Fixes
- **VOC Name Mapping**: Added comprehensive VOC name to GECKO code mapping
- **G-Prefix Normalization**: Fixed species code handling
- **Cytoscape JSON Format**: Fixed output format
- **SMILES Cleaning**: Fixed radical notation for browser compatibility
- **RDKit Rendering**: Replaced text boxes with molecular drawings
- **Restored 8 Missing Plots**: Van Krevelen, Top 10 species, SOA mass, etc.

---

## [2.0.0] - 2025-01-15

### Major Refactoring Release

#### Concurrency & Reliability
- **WorkspaceManager**: Isolated workspaces for each job
- **JobManager**: Thread-safe job state management with UUID-based IDs
- **SubprocessResult**: Structured subprocess execution

#### Chemical Informatics
- **MolecularGraph**: Graph-based molecular representation
- **GeckoGraphParser**: Proper GECKO notation parser
- **Known Species Database**: 100+ verified SMILES
- **RDKit Validation**: Optional SMILES validation

#### Scientific Rigor
- **PartitioningCalculator**: Dynamic gas-particle partitioning
- **VaporPressureReader**: Reads GECKO-A SAR vapor pressures
- **FortranOutputParser**: Smart fixed-width parsing

---

## [1.0.0] - Initial Release

### Features
- Basic GECKO-A integration
- Basic Box Model integration
- Web interface for job submission
- Simple result visualization

### Known Issues (Fixed in 2.0.0)
- Singleton input file causes race conditions
- Job ID collisions under load
- Regex-based SMILES fails for complex molecules
- Static C_OA assumption in partitioning
- MW-based vapor pressure estimation
- Fortran parsing errors on merged columns
