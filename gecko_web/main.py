"""
GECKO-A Web Interface - Main Application

Refactored for Production with:
- Proper concurrency handling (isolated workspaces)
- UUID-based job IDs
- Robust subprocess management
- Persistent state support

Author: Deeksha Sharma
"""

import os
import time
import shutil
import logging
import glob
import subprocess
import re
import json
import uuid
import tempfile
import threading
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, List, Any
from contextlib import contextmanager
from dataclasses import dataclass, field, asdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fastapi import FastAPI, Request, BackgroundTasks, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel

# Import internal modules
from gecko_web import mechanism_diagram
from gecko_web import postprocessing, reaction_tree
from gecko_web import pathway_visualizer
from gecko_web import mass_balance
from gecko_web import enhanced_visualizer
from gecko_web.chemdata import compound_database, voc_categories, reaction_data

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI(title="GECKO-A Interface", version="3.0.6")

# ==============================================================================
# Configuration
# ==============================================================================

BASE_DIR = Path(__file__).parent.absolute()
DATA_DIR = Path(os.getenv("DATA_DIR", "/data"))
if not DATA_DIR.exists():
    DATA_DIR = BASE_DIR.parent / "data"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES_DIR = BASE_DIR.parent / "samples"
GECKO_ROOT = Path(os.getenv("GECKO_SOURCE_DIR", "/app/gecko_source"))
BOXMODEL_ROOT = Path(os.getenv("BOXMODEL_SOURCE_DIR", "/app/boxmodel_source"))

# Workspace configuration
# Use a short path to avoid Fortran 100-char filename buffer overflow
# macOS temp dirs like /var/folders/6g/.../T/ are too long
WORKSPACE_BASE = Path("/tmp/gk")  # Short path for Fortran compatibility
if not WORKSPACE_BASE.exists():
    try:
        WORKSPACE_BASE.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        # Fall back to project directory if /tmp not writable
        WORKSPACE_BASE = BASE_DIR.parent / "workspaces"
        WORKSPACE_BASE.mkdir(parents=True, exist_ok=True)

# Job retention (seconds)
JOB_RETENTION_SECONDS = 86400  # 24 hours

# Mount static files and templates
app.mount("/static", StaticFiles(directory=str(BASE_DIR / "static")), name="static")
app.mount("/data", StaticFiles(directory=str(DATA_DIR)), name="data")
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))

# ==============================================================================
# Job State Management
# ==============================================================================

@dataclass
class JobState:
    """Represents the state of a job with full metadata."""
    id: str
    job_type: str
    voc_name: str
    status: str  # queued, running, completed, failed
    created_at: float
    updated_at: float
    logs: List[str] = field(default_factory=list)
    error: Optional[str] = None
    result_path: Optional[str] = None
    workspace_path: Optional[str] = None
    archived: bool = False
    archive_path: Optional[str] = None
    scenario_params: Optional[Dict] = None

    def to_dict(self) -> Dict:
        return asdict(self)


class JobManager:
    """
    Thread-safe job state manager with persistence support.
    """
    def __init__(self, storage_path: Optional[Path] = None):
        self._jobs: Dict[str, JobState] = {}
        self._lock = threading.RLock()
        self._storage_path = storage_path or (DATA_DIR / "job_state.json")
        self._load_state()

    def _load_state(self):
        """Load job state from disk on startup."""
        if self._storage_path.exists():
            try:
                with open(self._storage_path, 'r') as f:
                    data = json.load(f)
                    for job_id, job_data in data.items():
                        self._jobs[job_id] = JobState(**job_data)
                logger.info(f"Loaded {len(self._jobs)} jobs from persistent storage")
            except Exception as e:
                logger.error(f"Failed to load job state: {e}")

    def _save_state(self):
        """Persist job state to disk."""
        try:
            with open(self._storage_path, 'w') as f:
                data = {jid: j.to_dict() for jid, j in self._jobs.items()}
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save job state: {e}")

    def create_job(self, job_type: str, voc_name: str, scenario_params: Optional[Dict] = None) -> JobState:
        """Create a new job with a unique UUID-based ID."""
        with self._lock:
            job_id = f"{job_type}-{uuid.uuid4().hex[:12]}"
            now = time.time()
            job = JobState(
                id=job_id,
                job_type=job_type,
                voc_name=voc_name,
                status="queued",
                created_at=now,
                updated_at=now,
                scenario_params=scenario_params
            )
            self._jobs[job_id] = job
            self._save_state()
            return job

    def get_job(self, job_id: str) -> Optional[JobState]:
        with self._lock:
            return self._jobs.get(job_id)

    def update_job(self, job_id: str, **kwargs) -> Optional[JobState]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job:
                for key, value in kwargs.items():
                    if hasattr(job, key):
                        setattr(job, key, value)
                job.updated_at = time.time()
                self._save_state()
            return job

    def log_job(self, job_id: str, message: str):
        """Add a log message to a job."""
        with self._lock:
            job = self._jobs.get(job_id)
            if job:
                timestamp = datetime.now().strftime("%H:%M:%S")
                log_entry = f"[{timestamp}] {message}"
                job.logs.append(log_entry)
                job.updated_at = time.time()
                logger.info(f"[{job_id}] {message}")

    def list_jobs(self, include_old: bool = False) -> List[JobState]:
        """List all jobs, optionally filtering out old ones."""
        with self._lock:
            cutoff = time.time() - JOB_RETENTION_SECONDS if not include_old else 0
            return [j for j in self._jobs.values() if j.created_at >= cutoff]

    def delete_job(self, job_id: str) -> bool:
        with self._lock:
            if job_id in self._jobs:
                del self._jobs[job_id]
                self._save_state()
                return True
            return False

    def cleanup_old_jobs(self):
        """Remove jobs older than retention period."""
        with self._lock:
            cutoff = time.time() - JOB_RETENTION_SECONDS
            to_delete = [jid for jid, j in self._jobs.items()
                        if j.created_at < cutoff and j.status in ('completed', 'failed')]
            for jid in to_delete:
                del self._jobs[jid]
            if to_delete:
                self._save_state()
                logger.info(f"Cleaned up {len(to_delete)} old jobs")


# Global job manager instance
job_manager = JobManager()

# ==============================================================================
# Workspace Manager - Isolated Execution Environment
# ==============================================================================

class WorkspaceManager:
    """
    Manages isolated workspaces for GECKO-A execution.
    Each job gets its own workspace with:
    - Copy of GECKO executables and NML files
    - Isolated INPUT directory for cheminput.dat
    - Isolated OUT directory for results
    """

    def __init__(self, base_path: Path = WORKSPACE_BASE):
        self.base_path = base_path
        self.base_path.mkdir(parents=True, exist_ok=True)

    def create_workspace(self, job_id: str) -> Path:
        """Create an isolated workspace for a job."""
        # Use short directory name to avoid Fortran 100-char path limit
        # Extract just the UUID portion (last 12 chars of job_id)
        short_id = job_id.split('-')[-1] if '-' in job_id else job_id[:12]
        workspace = self.base_path / short_id

        if workspace.exists():
            shutil.rmtree(workspace)

        workspace.mkdir(parents=True)

        # Create directory structure
        (workspace / "INPUT").mkdir()
        (workspace / "RUN").mkdir()
        (workspace / "RUN" / "OUT").mkdir()

        return workspace

    def setup_gecko_workspace(self, workspace: Path, job_id: str) -> bool:
        """
        Copy GECKO executables and configuration to workspace.
        Returns True if successful, False otherwise.
        """
        try:
            # Copy essential files from GECKO source
            if GECKO_ROOT.exists():
                # Copy INPUT files (NML, databases)
                src_input = GECKO_ROOT / "INPUT"
                if src_input.exists():
                    for item in src_input.iterdir():
                        if item.is_file():
                            shutil.copy2(item, workspace / "INPUT")
                        elif item.is_dir():
                            shutil.copytree(item, workspace / "INPUT" / item.name, dirs_exist_ok=True)

                # Copy RUN scripts and executables
                src_run = GECKO_ROOT / "RUN"
                if src_run.exists():
                    for item in src_run.iterdir():
                        if item.is_file():
                            shutil.copy2(item, workspace / "RUN")
                        elif item.is_dir() and item.name != "OUT":
                            shutil.copytree(item, workspace / "RUN" / item.name, dirs_exist_ok=True)

                # Copy OBJ directory with cm executable (required by gecko.sh)
                src_obj = GECKO_ROOT / "OBJ"
                if src_obj.exists():
                    obj_dir = workspace / "OBJ"
                    obj_dir.mkdir(parents=True, exist_ok=True)
                    # Copy only the cm executable and essential files
                    cm_exe = src_obj / "cm"
                    if cm_exe.exists():
                        shutil.copy2(cm_exe, obj_dir / "cm")
                        # Make it executable
                        (obj_dir / "cm").chmod(0o755)
                    # Copy gecko.nml if it exists in OBJ
                    gecko_nml = src_obj / "gecko.nml"
                    if gecko_nml.exists():
                        shutil.copy2(gecko_nml, obj_dir / "gecko.nml")

                # Copy DATA directory (required for GECKO database files)
                src_data = GECKO_ROOT / "DATA"
                if src_data.exists():
                    shutil.copytree(src_data, workspace / "DATA", dirs_exist_ok=True)

                # Copy LIB directory if needed
                src_lib = GECKO_ROOT / "LIB"
                if src_lib.exists():
                    shutil.copytree(src_lib, workspace / "LIB", dirs_exist_ok=True)

                # Update paths in gecko.nml to point to workspace
                nml_path = workspace / "INPUT" / "gecko.nml"
                if nml_path.exists():
                    self._update_gecko_nml_paths(nml_path, workspace)

                return True
            else:
                logger.warning(f"GECKO_ROOT {GECKO_ROOT} does not exist")
                return False

        except Exception as e:
            logger.error(f"Failed to setup GECKO workspace: {e}")
            return False

    def _update_gecko_nml_paths(self, nml_path: Path, workspace: Path):
        """Update directory paths in gecko.nml to use workspace."""
        with open(nml_path, 'r') as f:
            content = f.read()

        workspace_str = str(workspace)
        if not workspace_str.endswith('/'):
            workspace_str += '/'

        out_path = str(workspace / "RUN" / "OUT")
        if not out_path.endswith('/'):
            out_path += '/'

        content = re.sub(r"dirgecko='.*?'", f"dirgecko='{workspace_str}'", content)
        content = re.sub(r"dirout='.*?'", f"dirout='{out_path}'", content)

        with open(nml_path, 'w') as f:
            f.write(content)

    def cleanup_workspace(self, workspace: Path, preserve: bool = False):
        """Remove workspace after job completion."""
        if not preserve and workspace.exists():
            try:
                shutil.rmtree(workspace)
                logger.debug(f"Cleaned up workspace: {workspace}")
            except Exception as e:
                logger.warning(f"Failed to cleanup workspace {workspace}: {e}")


workspace_manager = WorkspaceManager()

# ==============================================================================
# Subprocess Execution with Robust Error Handling
# ==============================================================================

@dataclass
class SubprocessResult:
    """Result of a subprocess execution."""
    success: bool
    return_code: int
    stdout: str
    stderr: str
    duration_seconds: float
    error_message: Optional[str] = None


def run_subprocess_robust(
    cmd: List[str],
    cwd: Path,
    description: str,
    timeout_seconds: int = 3600,
    job_id: Optional[str] = None
) -> SubprocessResult:
    """
    Execute a subprocess with robust error handling and logging.

    Args:
        cmd: Command and arguments to execute
        cwd: Working directory
        description: Human-readable description of the operation
        timeout_seconds: Maximum execution time
        job_id: Optional job ID for logging

    Returns:
        SubprocessResult with execution details
    """
    start_time = time.time()

    if job_id:
        job_manager.log_job(job_id, f"Executing: {' '.join(cmd)}")

    try:
        # Ensure script is executable
        if cmd[0].startswith('./'):
            script_path = cwd / cmd[0][2:]
            if script_path.exists():
                script_path.chmod(0o755)

        process = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ, 'LANG': 'C', 'LC_ALL': 'C'}
        )

        try:
            stdout, stderr = process.communicate(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()

            duration = time.time() - start_time
            error_msg = f"Process timed out after {timeout_seconds}s"

            if job_id:
                job_manager.log_job(job_id, f"TIMEOUT: {description}")

            return SubprocessResult(
                success=False,
                return_code=-1,
                stdout=stdout,
                stderr=stderr,
                duration_seconds=duration,
                error_message=error_msg
            )

        duration = time.time() - start_time
        success = process.returncode == 0

        # Log output (truncate if too long)
        if stdout and job_id:
            lines = stdout.strip().split('\n')
            log_lines = lines[-30:] if len(lines) > 30 else lines
            job_manager.log_job(job_id, f"STDOUT ({description}):\n" + '\n'.join(log_lines))

        if stderr and job_id:
            job_manager.log_job(job_id, f"STDERR ({description}):\n{stderr}")

        if success:
            if job_id:
                job_manager.log_job(job_id, f"Success: {description} ({duration:.1f}s)")
        else:
            error_msg = f"{description} failed with return code {process.returncode}"
            if job_id:
                job_manager.log_job(job_id, f"FAILED: {error_msg}")

        return SubprocessResult(
            success=success,
            return_code=process.returncode,
            stdout=stdout,
            stderr=stderr,
            duration_seconds=duration,
            error_message=None if success else f"Exit code {process.returncode}"
        )

    except FileNotFoundError as e:
        duration = time.time() - start_time
        error_msg = f"Command not found: {cmd[0]}"
        if job_id:
            job_manager.log_job(job_id, f"ERROR: {error_msg}")

        return SubprocessResult(
            success=False,
            return_code=-1,
            stdout="",
            stderr=str(e),
            duration_seconds=duration,
            error_message=error_msg
        )

    except Exception as e:
        duration = time.time() - start_time
        error_msg = f"Subprocess execution failed: {str(e)}"
        if job_id:
            job_manager.log_job(job_id, f"CRITICAL ERROR: {error_msg}")

        return SubprocessResult(
            success=False,
            return_code=-1,
            stdout="",
            stderr=str(e),
            duration_seconds=duration,
            error_message=error_msg
        )


# ==============================================================================
# VOC Mapping and Validation
# ==============================================================================

VOC_MAPPING = {
    # Biogenics
    "isoprene": "CH3Cd(=CdH2)CdH=CdH2",
    "alpha-pinene": "C12HCH2CH(C1(CH3)CH3)CH2CdH=Cd2CH3",
    "beta-pinene": "C12HCH2CH(C1(CH3)CH3)CH2CH2Cd2=CdH2",
    "myrcene": "CH3Cd(CH3)=CdHCH2CH2Cd(=CdH2)CdH=CdH2",
    "ocimene": "CH3Cd(CH3)=CdHCH2CdH=Cd(CH3)CdH=CdH2",
    "limonene": "C1H2CH2Cd(CH3)=CdHCH2C1HCd(CH3)=CdH2",

    # Alkanes
    "ethane": "CH3CH3",
    "propane": "CH3CH2CH3",
    "butane": "CH3CH2CH2CH3",
    "pentane": "CH3CH2CH2CH2CH3",
    "hexane": "CH3(CH2)4CH3",
    "heptane": "CH3(CH2)5CH3",
    "octane": "CH3(CH2)6CH3",
    "decane": "CH3(CH2)8CH3",
    "dodecane": "CH3(CH2)10CH3",
    "tetradecane": "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3",

    # Alkenes
    "ethene": "CH2=CH2",
    "propene": "CH3CH=CH2",
    "1-butene": "CH3CH2CH=CH2",
    "trans-2-butene": "CH3CdH=CdHCH3",

    # Aromatics
    "benzene": "c1HcHcHcHcHc1H",
    "toluene": "c1(CH3)cHcHcHcHc1H",
    "o-xylene": "c1(CH3)c(CH3)cHcHcHc1H",
    "m-xylene": "c1(CH3)cHc(CH3)cHcHc1H",
    "p-xylene": "c1(CH3)cHcHc(CH3)cHc1H",
    "ethylbenzene": "c1(CH2CH3)cHcHcHcHc1H"
}


def validate_voc(voc_name: str) -> bool:
    """Validate if the VOC name is plausible."""
    if not voc_name or len(voc_name) < 2:
        return False
    if voc_name.isdigit():
        return False
    # Check for dangerous characters
    if any(c in voc_name for c in [';', '|', '&', '$', '`', '\n', '\r']):
        return False
    return True


def get_gecko_input(voc_name: str) -> str:
    """
    Get GECKO-A input string for a VOC.

    This function checks multiple sources in order:
    1. Local VOC_MAPPING dictionary (legacy support)
    2. Compound database gecko_formula field
    3. Falls back to the VOC name itself

    Args:
        voc_name: Common name of the VOC (e.g., 'isoprene', 'cyclohexane')

    Returns:
        GECKO-A formatted input string
    """
    clean_name = voc_name.lower().strip()

    # Check local mapping first (for backward compatibility)
    if clean_name in VOC_MAPPING:
        return VOC_MAPPING[clean_name]

    # Check compound database
    try:
        from gecko_web.chemdata import get_gecko_formula
        gecko_formula = get_gecko_formula(clean_name)
        if gecko_formula:
            return gecko_formula
    except ImportError:
        pass

    # Try alternate name formats
    alt_names = [
        clean_name,
        clean_name.replace('-', '_'),
        clean_name.replace('_', '-'),
        clean_name.replace(' ', '_'),
        clean_name.replace(' ', '-'),
    ]

    try:
        from gecko_web.chemdata import get_compound
        for alt_name in alt_names:
            compound = get_compound(alt_name)
            if compound and compound.gecko_formula:
                return compound.gecko_formula
    except ImportError:
        pass

    # Return the name as-is (may work if it's already a valid GECKO formula)
    return voc_name


# ==============================================================================
# Pydantic Models
# ==============================================================================

class ScenarioParams(BaseModel):
    # Generator Params
    vapor_pressure_threshold: float = -13.0
    max_generations: int = 2

    # Box Model - Environment
    temperature_k: float = 298.0
    rh_percent: float = 50.0
    latitude_degrees: float = 45.0
    date_day: int = 21
    date_month: int = 6
    date_year: int = 1996

    # Box Model - Initial Conditions
    initial_o3_ppb: float = 40.0
    initial_nox_ppb: float = 10.0
    seed_aerosol_ug_m3: float = 10.0
    dilution_rate_s1: float = 0.0


class ExtendedGeckoOptions(BaseModel):
    """Extended GECKO-A Generator Options (Recommendation #8)"""

    # Vapor Pressure Method
    # Options: 'nannoolal', 'myrdal', 'simpol', 'evaporation'
    vapor_pressure_method: str = "nannoolal"

    # Critical vapor pressure threshold (log10 atm)
    # Species with Pvap < critvp are considered non-volatile
    critvp: float = -13.0

    # Maximum number of generations (oxidation steps)
    max_generations: int = 2

    # Reaction channel options
    enable_oh_reactions: bool = True
    enable_o3_reactions: bool = True
    enable_no3_reactions: bool = True
    enable_photolysis: bool = True
    enable_isomerization: bool = True

    # RO2 + RO2 permutation reactions
    enable_ro2_permutations: bool = True

    # PAN decomposition
    enable_pan_decomposition: bool = True

    # Nitrate formation
    # Options: 'carter', 'arey', 'atkinson'
    nitrate_yield_method: str = "carter"

    # Stereochemistry handling
    preserve_stereochemistry: bool = True

    # Ring-opening for cyclic peroxides
    enable_ring_opening: bool = True

    # SAR (Structure-Activity Relationship) options
    sar_oh_abstraction: str = "jenkin1997"  # or "atkinson2007"
    sar_oh_addition: str = "peeters2007"    # or "kwok1995"

    # Criegee intermediate chemistry
    enable_criegee_chemistry: bool = True
    sci_stabilization_fraction: float = 0.37  # Stabilized Criegee fraction

    # Autoxidation (H-shift reactions)
    enable_autoxidation: bool = False  # Computationally expensive
    autoxidation_rate_threshold: float = 0.01  # s^-1

    # Output options
    output_format: str = "gecko"  # 'gecko', 'kpp', 'mcm', 'facsimile'
    include_rate_coefficients: bool = True
    include_vapor_pressures: bool = True


class ExtendedBoxModelOptions(BaseModel):
    """Extended Box Model Options"""

    # Simulation time
    simulation_hours: float = 24.0
    output_interval_minutes: float = 5.0

    # Environment
    temperature_k: float = 298.0
    pressure_hpa: float = 1013.25
    rh_percent: float = 50.0

    # Location/Time
    latitude_degrees: float = 45.0
    longitude_degrees: float = 0.0
    date_day: int = 21
    date_month: int = 6
    date_year: int = 1996
    local_hour_start: float = 6.0  # Start time (local solar)

    # Initial Concentrations (ppb)
    initial_voc_ppb: float = 10.0
    initial_o3_ppb: float = 40.0
    initial_no_ppb: float = 5.0
    initial_no2_ppb: float = 5.0
    initial_hono_ppb: float = 0.1
    initial_h2o2_ppb: float = 1.0
    initial_hcho_ppb: float = 2.0
    initial_co_ppm: float = 0.1
    initial_ch4_ppm: float = 1.8

    # Aerosol parameters
    seed_aerosol_ug_m3: float = 10.0
    seed_aerosol_density_g_cm3: float = 1.4
    organic_aerosol_activity_coefficient: float = 1.0

    # Dilution/Deposition
    dilution_rate_s1: float = 0.0
    deposition_velocity_cm_s: float = 0.0

    # Chemistry
    enable_heterogeneous_reactions: bool = True
    enable_wall_loss: bool = False
    wall_loss_rate_s1: float = 0.0

    # Photolysis
    photolysis_scaling_factor: float = 1.0
    actinic_flux_scaling: float = 1.0


class VOCComparisonRequest(BaseModel):
    """Request for comparing multiple VOCs (Recommendation #14)"""
    voc_names: List[str]
    generator_options: Optional[ExtendedGeckoOptions] = None

    # Comparison metrics
    compare_mechanism_size: bool = True
    compare_product_distribution: bool = True
    compare_soa_yield: bool = True
    compare_radical_budget: bool = True


class MechanismReductionRequest(BaseModel):
    """Request for mechanism reduction (Recommendation #15)"""
    job_id: str  # Source job with full mechanism

    # Reduction method
    # Options: 'lumping', 'drgep', 'pfa', 'sensitivity'
    reduction_method: str = "drgep"

    # Target species to preserve
    target_species: List[str] = []  # Empty = auto-detect important species

    # Reduction parameters
    error_threshold: float = 0.1  # Max acceptable error (10%)
    min_species: int = 20  # Minimum species to keep
    preserve_radicals: bool = True
    preserve_soa_precursors: bool = True


class JobRequest(BaseModel):
    voc_name: str
    job_type: str  # 'generator', 'boxmodel', 'combined', 'comparison'
    scenario: Optional[ScenarioParams] = None
    extended_generator: Optional[ExtendedGeckoOptions] = None
    extended_boxmodel: Optional[ExtendedBoxModelOptions] = None


# ==============================================================================
# NML File Update Functions
# ==============================================================================

def update_gecko_nml(nml_path: Path, scenario: ScenarioParams):
    """Update gecko.nml with scenario parameters."""
    if not nml_path.exists():
        return

    with open(nml_path, 'r') as f:
        content = f.read()

    content = re.sub(r'critvp\s*=\s*[\d\.-]+', f'critvp = {scenario.vapor_pressure_threshold}', content)
    content = re.sub(r'maxgen\s*=\s*\d+', f'maxgen = {scenario.max_generations}', content)
    content = re.sub(r'TK\s*=\s*[\d\.]+', f'TK = {scenario.temperature_k}', content)

    with open(nml_path, 'w') as f:
        f.write(content)


def update_simu_nml(nml_path: Path, scenario: ScenarioParams):
    """Update simu.nml with scenario parameters."""
    if not nml_path.exists():
        return

    with open(nml_path, 'r') as f:
        content = f.read()

    content = re.sub(r'temp0\s*=\s*[\d\.]+', f'temp0 = {scenario.temperature_k}', content)
    content = re.sub(r'rh0\s*=\s*[\d\.]+', f'rh0 = {scenario.rh_percent}', content)
    content = re.sub(r'kdilu\s*=\s*[\d\.]+', f'kdilu = {scenario.dilution_rate_s1}', content)
    content = re.sub(r'sla\s*=\s*[\d\.]+', f'sla = {scenario.latitude_degrees}', content)
    content = re.sub(r'iday\s*=\s*\d+', f'iday = {scenario.date_day}', content)
    content = re.sub(r'imonth\s*=\s*\d+', f'imonth = {scenario.date_month}', content)
    content = re.sub(r'iyear\s*=\s*\d+', f'iyear = {scenario.date_year}', content)

    with open(nml_path, 'w') as f:
        f.write(content)


def update_concentrations_init(init_path: Path, scenario: ScenarioParams):
    """Update concentrations.init with initial conditions."""
    if not init_path.exists():
        return

    # Conversion: 1 ppb ~ 2.46e10 molecules/cm3 at 298K, 1 atm
    factor = 2.46e10
    o3_conc = scenario.initial_o3_ppb * factor
    no_conc = scenario.initial_nox_ppb * factor

    with open(init_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith('REAC GO3'):
            new_lines.append(f"REAC GO3          {o3_conc:.2E}\n")
        elif line.strip().startswith('REAC GNO '):
            new_lines.append(f"REAC GNO          {no_conc:.2E}\n")
        elif line.strip().startswith('REAC GNO2'):
            new_lines.append(f"!REAC GNO2        0.00E00\n")
        else:
            new_lines.append(line)

    with open(init_path, 'w') as f:
        f.writelines(new_lines)


# ==============================================================================
# GECKO-A Execution Functions
# ==============================================================================

def run_gecko_generator(job_id: str, voc_name: str, output_dir: Path,
                        scenario: Optional[ScenarioParams] = None) -> bool:
    """
    Run GECKO-A mechanism generator in an isolated workspace.

    Returns True on success, False on failure.
    """
    clean_name = voc_name.lower().strip()
    chem_input = get_gecko_input(voc_name)

    # Check library cache first
    library_path = DATA_DIR / "library" / clean_name
    if library_path.exists():
        job_manager.log_job(job_id, f"Found cached mechanism for {clean_name}")
        for item in library_path.iterdir():
            if item.is_file():
                shutil.copy2(item, output_dir)
        _create_reactions_file(output_dir)
        return True

    # Create isolated workspace
    job_manager.log_job(job_id, "Creating isolated workspace...")
    workspace = workspace_manager.create_workspace(job_id)

    try:
        # Setup GECKO in workspace
        if not workspace_manager.setup_gecko_workspace(workspace, job_id):
            job_manager.log_job(job_id, "WARNING: Could not setup GECKO workspace, trying direct execution")
            workspace = GECKO_ROOT
        else:
            job_manager.update_job(job_id, workspace_path=str(workspace))

        # Write input file to ISOLATED workspace
        input_file = workspace / "INPUT" / "cheminput.dat"
        with open(input_file, 'w') as f:
            f.write(f"{chem_input}\nEND\n")
        job_manager.log_job(job_id, f"Wrote input: {chem_input}")

        # Apply scenario parameters
        if scenario:
            nml_path = workspace / "INPUT" / "gecko.nml"
            update_gecko_nml(nml_path, scenario)
            job_manager.log_job(job_id, f"Applied scenario: VP={scenario.vapor_pressure_threshold}, MaxGen={scenario.max_generations}")

        # Execute GECKO-A
        run_dir = workspace / "RUN"
        result = run_subprocess_robust(
            ["./gecko.sh"],
            run_dir,
            "GECKO-A Mechanism Generator",
            timeout_seconds=3600,
            job_id=job_id
        )

        if not result.success:
            raise RuntimeError(f"GECKO-A failed: {result.error_message}")

        # Copy results from workspace OUT to job output
        out_dir = workspace / "RUN" / "OUT"
        if out_dir.exists():
            for item in out_dir.iterdir():
                if item.is_file():
                    shutil.copy2(item, output_dir)
            job_manager.log_job(job_id, f"Copied {len(list(out_dir.iterdir()))} result files")

        _create_reactions_file(output_dir)
        return True

    except Exception as e:
        job_manager.log_job(job_id, f"GECKO-A execution failed: {str(e)}")
        raise
    finally:
        # Cleanup workspace (preserve on failure for debugging in dev)
        if os.getenv("PRESERVE_WORKSPACES", "").lower() != "true":
            workspace_manager.cleanup_workspace(workspace)


def _create_reactions_file(output_dir: Path):
    """Create consolidated reactions.txt file from mechanism outputs."""
    reaction_file = output_dir / "reactions.txt"

    potential_files = (
        list(output_dir.glob("*.mec")) +
        list(output_dir.glob("*.k")) +
        list(output_dir.glob("*.mech")) +
        list(output_dir.glob("*dictionary*"))
    )

    if potential_files:
        with open(reaction_file, 'w') as outfile:
            for pf in potential_files:
                outfile.write(f"--- {pf.name} ---\n")
                with open(pf, 'r') as infile:
                    outfile.write(infile.read())
                outfile.write("\n\n")
    else:
        with open(reaction_file, 'w') as f:
            f.write("GECKO-A ran, but no mechanism file was found.\n")
            f.write("Files generated:\n")
            f.write("\n".join(item.name for item in output_dir.iterdir()))


def run_box_model(job_id: str, voc_name: str, output_dir: Path,
                  scenario: Optional[ScenarioParams] = None) -> bool:
    """
    Run BOXMODEL4GECKO simulation.

    Returns True on success, raises exception on failure.
    """
    if not BOXMODEL_ROOT.exists():
        raise EnvironmentError(
            f"BOXMODEL4GECKO not found at {BOXMODEL_ROOT}. "
            "Please run in Docker with proper environment."
        )

    prepare_script = BOXMODEL_ROOT / "prepare_simu.sh"
    if not prepare_script.exists():
        raise EnvironmentError(f"Box Model prepare_simu.sh not found at {prepare_script}")

    mech_name = voc_name.lower().strip()
    mech_upper = mech_name.upper()

    # Clean up previous simulation
    simu_dir = BOXMODEL_ROOT / "SIMU" / mech_upper
    if simu_dir.exists():
        job_manager.log_job(job_id, f"Removing existing simulation: {simu_dir}")
        shutil.rmtree(simu_dir)

    # Prepare simulation
    job_manager.log_job(job_id, f"Preparing Box Model for {voc_name}...")
    result = run_subprocess_robust(
        ["./prepare_simu.sh", "--import_mech", mech_name, "--geckooutdir", str(output_dir)],
        BOXMODEL_ROOT,
        "Box Model Preparation",
        job_id=job_id
    )

    if not result.success:
        raise RuntimeError(f"Box Model preparation failed: {result.error_message}")

    # Apply scenario parameters
    simu_dir = BOXMODEL_ROOT / "SIMU" / mech_upper
    if scenario and simu_dir.exists():
        job_manager.log_job(job_id, "Applying scenario parameters to simulation...")
        update_simu_nml(simu_dir / "simu.nml", scenario)
        update_concentrations_init(simu_dir / "concentrations.init", scenario)

    # Run simulation
    start_script = simu_dir / "start.sh"
    if not start_script.exists():
        raise RuntimeError(f"Simulation directory not created properly: {simu_dir}")

    job_manager.log_job(job_id, f"Running Box Model simulation...")
    result = run_subprocess_robust(
        ["./start.sh"],
        simu_dir,
        "Box Model Execution",
        timeout_seconds=7200,
        job_id=job_id
    )

    if not result.success:
        raise RuntimeError(f"Box Model execution failed: {result.error_message}")

    # Copy results
    resu_dir = simu_dir / "RESU"
    if resu_dir.exists():
        for item in resu_dir.iterdir():
            if item.is_file():
                shutil.copy2(item, output_dir)

    # Copy mechanism file from CHEMDAT
    mech_file = BOXMODEL_ROOT / "CHEMDAT" / f"{mech_name}.mech"
    if mech_file.exists():
        shutil.copy2(mech_file, output_dir / f"{mech_name}.mech")

    return True


def run_box_model_with_mechanism_dir(
    job_id: str,
    voc_name: str,
    mechanism_dir: Path,
    output_dir: Path,
    scenario: Optional[ScenarioParams] = None
) -> bool:
    """
    Run BOXMODEL4GECKO simulation with explicit mechanism and output directories.

    This function is designed for the combined workflow where mechanism files
    are stored in a separate directory from the simulation output.

    Args:
        job_id: Job ID for logging
        voc_name: VOC name
        mechanism_dir: Directory containing mechanism files (dictionary.out, *.mec, etc.)
        output_dir: Directory to store simulation results
        scenario: Optional simulation parameters

    Returns:
        True on success, raises exception on failure.
    """
    if not BOXMODEL_ROOT.exists():
        raise EnvironmentError(
            f"BOXMODEL4GECKO not found at {BOXMODEL_ROOT}. "
            "Please run in Docker with proper environment."
        )

    prepare_script = BOXMODEL_ROOT / "prepare_simu.sh"
    if not prepare_script.exists():
        raise EnvironmentError(f"Box Model prepare_simu.sh not found at {prepare_script}")

    # Verify mechanism files exist
    mechanism_dir = Path(mechanism_dir)
    if not mechanism_dir.exists():
        raise FileNotFoundError(f"Mechanism directory not found: {mechanism_dir}")

    # Check for required mechanism files
    dict_file = mechanism_dir / "dictionary.out"
    if not dict_file.exists():
        # List available files for debugging
        available = [f.name for f in mechanism_dir.iterdir() if f.is_file()]
        raise FileNotFoundError(
            f"dictionary.out not found in {mechanism_dir}. "
            f"Available files: {available[:10]}"
        )

    job_manager.log_job(job_id, f"Using mechanism from: {mechanism_dir}")

    mech_name = voc_name.lower().strip()
    mech_upper = mech_name.upper()

    # Clean up previous simulation
    simu_dir = BOXMODEL_ROOT / "SIMU" / mech_upper
    if simu_dir.exists():
        job_manager.log_job(job_id, f"Removing existing simulation: {simu_dir}")
        shutil.rmtree(simu_dir)

    # Prepare simulation - use mechanism_dir for geckooutdir
    job_manager.log_job(job_id, f"Preparing Box Model for {voc_name}...")
    result = run_subprocess_robust(
        ["./prepare_simu.sh", "--import_mech", mech_name, "--geckooutdir", str(mechanism_dir)],
        BOXMODEL_ROOT,
        "Box Model Preparation",
        job_id=job_id
    )

    if not result.success:
        # Log more details about the failure
        job_manager.log_job(job_id, f"prepare_simu.sh failed. Checking mechanism files...")
        mec_files = list(mechanism_dir.glob("*.mec"))
        job_manager.log_job(job_id, f"Found {len(mec_files)} .mec files: {[f.name for f in mec_files[:5]]}")
        raise RuntimeError(f"Box Model preparation failed: {result.error_message}")

    # Apply scenario parameters
    simu_dir = BOXMODEL_ROOT / "SIMU" / mech_upper
    if not simu_dir.exists():
        raise RuntimeError(
            f"Simulation directory not created: {simu_dir}. "
            "The prepare_simu.sh script may have failed silently."
        )

    if scenario:
        job_manager.log_job(job_id, "Applying scenario parameters to simulation...")
        update_simu_nml(simu_dir / "simu.nml", scenario)
        update_concentrations_init(simu_dir / "concentrations.init", scenario)

    # Run simulation
    start_script = simu_dir / "start.sh"
    if not start_script.exists():
        raise RuntimeError(f"start.sh not found in simulation directory: {simu_dir}")

    job_manager.log_job(job_id, f"Running Box Model simulation...")
    result = run_subprocess_robust(
        ["./start.sh"],
        simu_dir,
        "Box Model Execution",
        timeout_seconds=7200,
        job_id=job_id
    )

    if not result.success:
        raise RuntimeError(f"Box Model execution failed: {result.error_message}")

    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy results to output directory
    resu_dir = simu_dir / "RESU"
    if resu_dir.exists():
        result_files = list(resu_dir.iterdir())
        job_manager.log_job(job_id, f"Copying {len(result_files)} result files to output...")
        for item in result_files:
            if item.is_file():
                shutil.copy2(item, output_dir)
    else:
        job_manager.log_job(job_id, "WARNING: No RESU directory found after simulation")

    # Copy mechanism file from CHEMDAT
    mech_file = BOXMODEL_ROOT / "CHEMDAT" / f"{mech_name}.mech"
    if mech_file.exists():
        shutil.copy2(mech_file, output_dir / f"{mech_name}.mech")

    job_manager.log_job(job_id, "Box model simulation completed successfully")
    return True


# ==============================================================================
# Main Job Processing
# ==============================================================================

def run_process(job_id: str, request: JobRequest):
    """Main job processing function."""
    voc_name = request.voc_name
    job_type = request.job_type
    scenario = request.scenario

    job_manager.log_job(job_id, f"Starting {job_type} job for {voc_name}")
    job_manager.update_job(job_id, status='running')

    try:
        if not validate_voc(voc_name):
            raise ValueError(f"Invalid VOC name: '{voc_name}'")

        output_dir = DATA_DIR / "output" / job_id
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run Generator
        if job_type == 'generator' or not (output_dir / "dictionary.out").exists():
            job_manager.log_job(job_id, "Running GECKO-A generator...")
            run_gecko_generator(job_id, voc_name, output_dir, scenario)

        if job_type == 'generator':
            # Generate diagrams for generator job
            job_manager.log_job(job_id, "Generating mechanism diagrams...")
            try:
                # Always use robust parser from reaction_tree.py to get functional groups and SMILES
                tree_data = reaction_tree.parse_reaction_tree(str(output_dir), voc_name, max_depth=3, max_nodes=50)
                json_path = output_dir / "reaction_tree.json"
                with open(json_path, 'w') as f:
                    json.dump(tree_data, f)
                
                # Generate enhanced diagrams
                visualizer = enhanced_visualizer.EnhancedPathwayVisualizer()
                visualizer.from_tree_data(tree_data)
                
                visualizer.generate_diagram(
                    output_path=str(output_dir / "pathway_diagram"),
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    layout='hierarchical',
                    show_branching_ratios=True,
                    color_by_reaction_type=True,
                    format='png'
                )
                visualizer.generate_diagram(
                    output_path=str(output_dir / "pathway_diagram"),
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    layout='hierarchical',
                    show_branching_ratios=True,
                    color_by_reaction_type=True,
                    format='svg'
                )
                job_manager.log_job(job_id, "Diagrams generated successfully.")
            except Exception as e:
                job_manager.log_job(job_id, f"Diagram generation failed: {e}")

            job_manager.log_job(job_id, "Generator job completed.")
            job_manager.update_job(
                job_id,
                status='completed',
                result_path=f"/data/output/{job_id}"
            )
            return

        # Run Box Model
        job_manager.log_job(job_id, "Running Box Model simulation...")
        run_box_model(job_id, voc_name, output_dir, scenario)

        # Post-processing
        job_manager.log_job(job_id, "Running post-processing...")
        try:
            postprocessing.run_postprocessing(str(output_dir), voc_name)
        except Exception as e:
            job_manager.log_job(job_id, f"Post-processing warning: {e}")

        # Recreate reactions.txt with all files
        _create_reactions_file(output_dir)

        # Generate mechanism diagrams
        job_manager.log_job(job_id, "Generating mechanism diagrams...")
        try:
            mechanism_diagram.generate_all_outputs(
                output_dir=str(output_dir),
                voc_name=voc_name,
                generate_static=True,
                generate_kpp=True,
                generate_mcm=True,
                generate_facsimile=True,
                max_depth=3,
                min_branching=0.05,
                max_nodes=25,
                max_children=4
            )
        except Exception as e:
            job_manager.log_job(job_id, f"Diagram generation warning: {e}")
            try:
                tree_data = reaction_tree.parse_reaction_tree(str(output_dir), voc_name, max_depth=3, max_nodes=30)
                with open(output_dir / "reaction_tree.json", 'w') as f:
                    json.dump(tree_data, f)
            except Exception as e2:
                job_manager.log_job(job_id, f"Reaction tree fallback failed: {e2}")

        # Generate publication-quality pathway diagram using enhanced visualizer
        job_manager.log_job(job_id, "Generating publication-quality pathway diagram...")
        try:
            tree_file = output_dir / "reaction_tree.json"
            if tree_file.exists():
                with open(tree_file) as f:
                    tree_data = json.load(f)

                # Generate diagrams using Enhanced Visualizer
                visualizer = enhanced_visualizer.EnhancedPathwayVisualizer()
                visualizer.from_tree_data(tree_data)

                # Generate PNG
                visualizer.generate_diagram(
                    output_path=str(output_dir / "pathway_diagram"),
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    layout='hierarchical',
                    show_branching_ratios=True,
                    color_by_reaction_type=True,
                    format='png'
                )
                
                # Generate SVG
                visualizer.generate_diagram(
                    output_path=str(output_dir / "pathway_diagram"),
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    layout='hierarchical',
                    show_branching_ratios=True,
                    color_by_reaction_type=True,
                    format='svg'
                )
                
                job_manager.log_job(job_id, "Enhanced pathway diagram generated successfully.")
        except Exception as e:
            job_manager.log_job(job_id, f"Enhanced visualizer failed ({e}), falling back to legacy...")
            try:
                 pathway_visualizer.generate_pathway_diagram(
                    tree_data=tree_data,
                    output_path=str(output_dir / "pathway_diagram"),
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    format='png',
                    mol_size=(180, 140)
                )
            except Exception as e2:
                 job_manager.log_job(job_id, f"Pathway diagram generation completely failed: {e2}")

        job_manager.log_job(job_id, f"[{job_type.upper()}] Completed successfully.")
        job_manager.update_job(
            job_id,
            status='completed',
            result_path=f"/data/output/{job_id}"
        )

    except Exception as e:
        job_manager.log_job(job_id, f"Job failed: {str(e)}")
        job_manager.update_job(job_id, status='failed', error=str(e))


# ==============================================================================
# API Endpoints
# ==============================================================================

@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/api/environment")
async def check_environment():
    """Check if GECKO-A and Box Model are available."""
    gecko_available = (GECKO_ROOT / "RUN" / "gecko.sh").exists()
    boxmodel_available = BOXMODEL_ROOT.exists() and (BOXMODEL_ROOT / "prepare_simu.sh").exists()

    errors = []
    if not gecko_available:
        errors.append(f"GECKO-A not found at {GECKO_ROOT}")
    if not boxmodel_available:
        errors.append(f"BOXMODEL not found at {BOXMODEL_ROOT}")

    return {
        "gecko_available": gecko_available,
        "boxmodel_available": boxmodel_available,
        "fully_functional": gecko_available and boxmodel_available,
        "errors": errors,
        "message": (
            "Environment fully configured."
            if (gecko_available and boxmodel_available)
            else "Running in limited mode."
        )
    }


@app.post("/api/jobs")
async def create_job(job_req: JobRequest, background_tasks: BackgroundTasks):
    """Create and queue a new job."""
    scenario_dict = job_req.scenario.dict() if job_req.scenario else None
    job = job_manager.create_job(job_req.job_type, job_req.voc_name, scenario_dict)

    background_tasks.add_task(run_process, job.id, job_req)

    return {"job_id": job.id, "status": "queued"}


@app.get("/api/jobs/{job_id}")
async def get_job_status(job_id: str):
    """Get status of a specific job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job.to_dict()


@app.get("/api/jobs/{job_id}/results")
async def get_job_results(job_id: str):
    """Get results of a completed job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != 'completed':
        raise HTTPException(status_code=400, detail="Job not completed")

    output_dir = DATA_DIR / "output" / job_id

    results = {
        "reactions": "",
        "csv_data": [],
        "plots": [],
        "reaction_tree": None,
        "mechanism_exports": [],
        "mechanism_summary": None,
        "static_diagram": None,
        "smiles_data": []
    }

    # Read reactions
    reactions_file = output_dir / "reactions.txt"
    if reactions_file.exists():
        results["reactions"] = reactions_file.read_text()

    # Read reaction tree
    tree_file = output_dir / "reaction_tree.json"
    if tree_file.exists():
        with open(tree_file) as f:
            tree_data = json.load(f)
            results["reaction_tree"] = tree_data

            if tree_data and "nodes" in tree_data:
                results["smiles_data"] = [
                    {
                        "code": node.get("display_code", node["id"]),
                        "smiles": node["smiles"],
                        "formula": node.get("label", ""),
                        "molecular_weight": node.get("molecular_weight", 0),
                        "functional_groups": node.get("functional_groups", [])
                    }
                    for node in tree_data["nodes"]
                    if node.get("smiles") and node["smiles"] != node["id"]
                ]

    # Mechanism summary
    summary_file = output_dir / "mechanism_summary.json"
    if summary_file.exists():
        with open(summary_file) as f:
            results["mechanism_summary"] = json.load(f)

    # Pre-calculate paths for diagrams
    pathway_png = output_dir / "pathway_diagram.png"
    pathway_svg = output_dir / "pathway_diagram.svg"
    diagram_png = output_dir / "mechanism_diagram.png"
    diagram_svg = output_dir / "mechanism_diagram.svg"

    # Mechanism exports
    export_map = {'.kpp': 'KPP', '.mcm': 'MCM', '.fac': 'FACSIMILE'}
    for item in output_dir.iterdir():
        ext = item.suffix.lower()
        if ext in export_map:
            results["mechanism_exports"].append({
                "format": export_map[ext],
                "filename": item.name,
                "url": f"/data/output/{job_id}/{item.name}"
            })

    # CSV data
    csv_file = output_dir / "aerosol_data.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
        results["csv_data"] = df.head(100).to_dict(orient="records")

    # Result filtering: exclude diagrams from the general plots list
    # Because they are shown in dedicated sections
    excluded_plots = {'pathway_diagram.png', 'mechanism_diagram.png'}

    # Plots
    for item in sorted(output_dir.iterdir()):
        if item.suffix == ".png" and item.name not in excluded_plots:
            results["plots"].append({
                "title": item.stem.replace("_", " ").title(),
                "url": f"/data/output/{job_id}/{item.name}"
            })
    
    # NEW: Also include the enhanced diagram explicitly for the frontend
    # This allows the frontend to show download buttons for the main diagram
    if pathway_png.exists():
        results["pathway_diagram"] = {
            "png": f"/data/output/{job_id}/pathway_diagram.png",
            "svg": f"/data/output/{job_id}/pathway_diagram.svg" if pathway_svg.exists() else None
        }

    # Static diagram (try enhanced first, then legacy)
    # This maintains the legacy "Mechanism Diagram" view (can be same content)
    if diagram_png.exists():
        results["static_diagram"] = {
            "png": f"/data/output/{job_id}/mechanism_diagram.png",
            "svg": f"/data/output/{job_id}/mechanism_diagram.svg" if diagram_svg.exists() else None
        }
    elif pathway_png.exists():
        # Fallback: Populate static_diagram slot with pathway_diagram if legacy doesn't exist
        results["static_diagram"] = {
            "png": f"/data/output/{job_id}/pathway_diagram.png",
            "svg": f"/data/output/{job_id}/pathway_diagram.svg" if pathway_svg.exists() else None
        }

    return results


@app.post("/api/jobs/{job_id}/archive")
async def archive_job(job_id: str):
    """Archive a completed job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != 'completed':
        raise HTTPException(status_code=400, detail="Job not completed")

    source_dir = DATA_DIR / "output" / job_id
    if not source_dir.exists():
        raise HTTPException(status_code=404, detail="Job data not found")

    safe_voc = re.sub(r'[^a-zA-Z0-9_-]', '_', job.voc_name)
    date_str = datetime.fromtimestamp(job.created_at).strftime('%Y-%m-%d')
    archive_dir = DATA_DIR / "archives" / safe_voc / date_str / job_id
    archive_dir.mkdir(parents=True, exist_ok=True)

    try:
        for item in source_dir.iterdir():
            if item.is_file():
                shutil.copy2(item, archive_dir)
            elif item.is_dir():
                shutil.copytree(item, archive_dir / item.name, dirs_exist_ok=True)

        job_manager.update_job(job_id, archived=True, archive_path=str(archive_dir))
        return {"message": "Job archived successfully", "path": str(archive_dir)}

    except Exception as e:
        logger.error(f"Failed to archive job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to archive: {str(e)}")


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if not job.archived:
        output_dir = DATA_DIR / "output" / job_id
        if output_dir.exists():
            try:
                shutil.rmtree(output_dir)
            except Exception as e:
                logger.error(f"Failed to delete output for {job_id}: {e}")

    job_manager.delete_job(job_id)
    return {"message": "Job deleted successfully"}


@app.get("/api/jobs")
async def list_jobs():
    """List all jobs."""
    return [j.to_dict() for j in job_manager.list_jobs()]


# ==============================================================================
# Extended API Endpoints - Compound Database
# ==============================================================================

@app.get("/api/compounds")
async def list_compounds():
    """List all available VOC compounds with their properties."""
    compounds = []
    for name, info in compound_database.COMPOUND_DATABASE.items():
        compounds.append({
            "name": name,
            "smiles": info.smiles,
            "gecko_formula": info.gecko_formula,
            "category": info.category,
            "subcategory": info.subcategory,
            "molecular_weight": info.molecular_weight,
            "molecular_formula": info.molecular_formula,
            "vapor_pressure_pa": info.vapor_pressure_298k_pa,
            "oh_rate_constant": info.koh_298k,
            "o3_rate_constant": info.ko3_298k,
            "no3_rate_constant": info.kno3_298k,
            "boiling_point_k": info.boiling_point_k,
            "cas": info.cas
        })
    return {"compounds": compounds, "count": len(compounds)}


@app.get("/api/compounds/categories")
async def list_compound_categories():
    """List all VOC categories with their compounds."""
    categories = {}
    for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
        categories[cat_id] = {
            "name": cat_info.name,
            "description": cat_info.description,
            "subcategories": {
                sub_id: {
                    "name": sub_info.name,
                    "description": sub_info.description,
                    "compounds": sub_info.compounds
                }
                for sub_id, sub_info in cat_info.subcategories.items()
            }
        }
    return {"categories": categories}


@app.get("/api/compounds/dropdown")
async def get_compounds_for_dropdown():
    """Fast endpoint optimized for populating the compound dropdown."""
    # Build a simple structure: { category_name: [compound_names] }
    dropdown_data = {}
    total_count = 0

    for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
        compounds_in_category = []
        for sub_id, sub_info in cat_info.subcategories.items():
            compounds_in_category.extend(sub_info.compounds)

        # Remove duplicates and sort
        compounds_in_category = sorted(set(compounds_in_category))
        if compounds_in_category:
            dropdown_data[cat_info.name] = compounds_in_category
            total_count += len(compounds_in_category)

    return {
        "categories": dropdown_data,
        "total_compounds": total_count
    }


@app.get("/api/compounds/{compound_name}")
async def get_compound_detail(compound_name: str):
    """Get detailed information about a specific compound."""
    compound = compound_database.get_compound(compound_name)
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound '{compound_name}' not found")

    return {
        "name": compound.name,
        "smiles": compound.smiles,
        "gecko_formula": compound.gecko_formula,
        "molecular_formula": compound.molecular_formula,
        "category": compound.category,
        "subcategory": compound.subcategory,
        "molecular_weight": compound.molecular_weight,
        "vapor_pressure_pa": compound.vapor_pressure_298k_pa,
        "boiling_point_k": compound.boiling_point_k,
        "henrys_law_constant": compound.henrys_law_mol_m3_pa,
        "oh_rate_constant": compound.koh_298k,
        "o3_rate_constant": compound.ko3_298k,
        "no3_rate_constant": compound.kno3_298k,
        "atmospheric_lifetime_hours": compound_database.get_atmospheric_lifetime(compound_name),
        "cas": compound.cas,
        "inchi": compound.inchi,
        "aliases": compound.aliases,
        "notes": compound.notes,
        "source": compound.source
    }


@app.get("/api/compounds/search/{query}")
async def search_compounds(query: str):
    """Search compounds by name, category, or property."""
    results = compound_database.search_compounds(query)
    return {"query": query, "results": results, "count": len(results)}


@app.get("/api/compounds/category/{category}")
async def get_compounds_by_category(category: str):
    """Get all compounds in a specific category."""
    compounds = voc_categories.get_category_compounds(category)
    if not compounds:
        raise HTTPException(status_code=404, detail=f"Category '{category}' not found")
    return {"category": category, "compounds": compounds, "count": len(compounds)}


# ==============================================================================
# Extended API Endpoints - Reaction Kinetics
# ==============================================================================

@app.get("/api/kinetics/{compound_name}")
async def get_reaction_kinetics(compound_name: str, temperature_k: float = 298.0):
    """Get reaction rate constants for a compound at specified temperature."""
    rates = {}

    # Try to get kinetics from database using correct method names
    oxidant_map = {
        "OH": reaction_data.OxidantType.OH,
        "O3": reaction_data.OxidantType.O3,
        "NO3": reaction_data.OxidantType.NO3
    }

    for oxidant_name, oxidant_type in oxidant_map.items():
        kinetics = reaction_data.reaction_database.get_reaction(compound_name, oxidant_type)
        if kinetics:
            rate = kinetics.rate_params.calculate_k(temperature_k)
            rates[oxidant_name] = {
                "rate_constant": rate,
                "temperature_k": temperature_k,
                "units": "cm3/molecule/s",
                "reference": kinetics.reference or kinetics.rate_params.reference
            }

    if not rates:
        # Fall back to compound database
        compound = compound_database.get_compound(compound_name)
        if compound:
            if compound.koh_298k and compound.koh_298k > 0:
                rates["OH"] = {
                    "rate_constant": compound.koh_298k,
                    "temperature_k": 298.0,
                    "units": "cm3/molecule/s",
                    "reference": "IUPAC/JPL recommended"
                }
            if compound.ko3_298k and compound.ko3_298k > 0:
                rates["O3"] = {
                    "rate_constant": compound.ko3_298k,
                    "temperature_k": 298.0,
                    "units": "cm3/molecule/s",
                    "reference": "IUPAC/JPL recommended"
                }
            if compound.kno3_298k and compound.kno3_298k > 0:
                rates["NO3"] = {
                    "rate_constant": compound.kno3_298k,
                    "temperature_k": 298.0,
                    "units": "cm3/molecule/s",
                    "reference": "IUPAC/JPL recommended"
                }

    if not rates:
        raise HTTPException(status_code=404, detail=f"No kinetics data for '{compound_name}'")

    return {
        "compound": compound_name,
        "temperature_k": temperature_k,
        "rate_constants": rates
    }


@app.get("/api/kinetics/{compound_name}/lifetime")
async def get_atmospheric_lifetime(
    compound_name: str,
    oh_concentration: float = 1e6,
    o3_concentration: float = 7e11,
    no3_concentration: float = 5e8,
    temperature_k: float = 298.0
):
    """Calculate atmospheric lifetime of a compound."""
    # Use the module-level function, not a method
    lifetimes = reaction_data.calculate_atmospheric_lifetime(
        compound_name,
        T=temperature_k,
        oh_conc=oh_concentration,
        o3_conc=o3_concentration,
        no3_conc=no3_concentration,
        db=reaction_data.reaction_database
    )

    if not lifetimes or all(v == float('inf') for v in lifetimes.values()):
        raise HTTPException(status_code=404, detail=f"No kinetics data for '{compound_name}'")

    return {
        "compound": compound_name,
        "conditions": {
            "OH_concentration_molec_cm3": oh_concentration,
            "O3_concentration_molec_cm3": o3_concentration,
            "NO3_concentration_molec_cm3": no3_concentration,
            "temperature_k": temperature_k
        },
        "lifetimes": lifetimes
    }


# ==============================================================================
# Extended API Endpoints - Mass Balance Verification
# ==============================================================================

@app.get("/api/jobs/{job_id}/mass-balance")
async def get_mass_balance(job_id: str):
    """Get mass balance verification results for a completed job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != 'completed':
        raise HTTPException(status_code=400, detail="Job not completed")

    output_dir = DATA_DIR / "output" / job_id

    # Check for existing mass balance report
    mb_file = output_dir / "mass_balance_report.json"
    if mb_file.exists():
        with open(mb_file) as f:
            return json.load(f)

    # Generate mass balance verification
    try:
        checker = mass_balance.MassBalanceChecker()
        dictionary_file = output_dir / "dictionary.out"
        mechanism_file = None

        # Find mechanism file
        for ext in ['.mec', '.mech', '.k']:
            candidates = list(output_dir.glob(f"*{ext}"))
            if candidates:
                mechanism_file = candidates[0]
                break

        if not dictionary_file.exists():
            raise HTTPException(status_code=404, detail="Dictionary file not found")

        # Parse and verify
        checker.parse_dictionary(dictionary_file)
        if mechanism_file:
            checker.parse_mechanism(mechanism_file)

        result = checker.verify_all_reactions()

        # Save for future requests
        with open(mb_file, 'w') as f:
            json.dump(result, f, indent=2)

        return result

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Mass balance verification failed: {str(e)}")


# ==============================================================================
# Extended API Endpoints - VOC Comparison (Recommendation #14)
# ==============================================================================

@app.post("/api/workflow/comparison")
async def create_comparison_workflow(
    request: VOCComparisonRequest,
    background_tasks: BackgroundTasks
):
    """Create a workflow to compare multiple VOCs."""
    if len(request.voc_names) < 2:
        raise HTTPException(status_code=400, detail="At least 2 VOCs required for comparison")
    if len(request.voc_names) > 10:
        raise HTTPException(status_code=400, detail="Maximum 10 VOCs for comparison")

    job = job_manager.create_job("comparison", ",".join(request.voc_names), {
        "voc_names": request.voc_names,
        "generator_options": request.generator_options.dict() if request.generator_options else None,
        "comparison_metrics": {
            "mechanism_size": request.compare_mechanism_size,
            "product_distribution": request.compare_product_distribution,
            "soa_yield": request.compare_soa_yield,
            "radical_budget": request.compare_radical_budget
        }
    })

    background_tasks.add_task(run_comparison_workflow, job.id, request)
    return {"job_id": job.id, "status": "queued", "workflow": "comparison", "vocs": request.voc_names}


async def run_comparison_workflow(job_id: str, request: VOCComparisonRequest):
    """Execute VOC comparison workflow."""
    job_manager.log_job(job_id, f"Starting comparison workflow for {len(request.voc_names)} VOCs")
    job_manager.update_job(job_id, status='running')

    try:
        output_dir = DATA_DIR / "output" / job_id
        output_dir.mkdir(parents=True, exist_ok=True)

        comparison_results = {
            "vocs": request.voc_names,
            "individual_results": {},
            "comparison": {}
        }

        # Process each VOC
        for i, voc_name in enumerate(request.voc_names):
            job_manager.log_job(job_id, f"Processing VOC {i+1}/{len(request.voc_names)}: {voc_name}")

            voc_dir = output_dir / voc_name.lower().replace(" ", "_")
            voc_dir.mkdir(exist_ok=True)

            # Run generator for this VOC
            try:
                scenario = ScenarioParams()
                if request.generator_options:
                    scenario.vapor_pressure_threshold = request.generator_options.critvp
                    scenario.max_generations = request.generator_options.max_generations

                run_gecko_generator(job_id, voc_name, voc_dir, scenario)

                # Collect mechanism statistics
                stats = _collect_mechanism_stats(voc_dir, voc_name)
                comparison_results["individual_results"][voc_name] = stats

            except Exception as e:
                job_manager.log_job(job_id, f"Failed to process {voc_name}: {e}")
                comparison_results["individual_results"][voc_name] = {"error": str(e)}

        # Generate comparison metrics
        if request.compare_mechanism_size:
            comparison_results["comparison"]["mechanism_size"] = _compare_mechanism_sizes(
                comparison_results["individual_results"]
            )

        if request.compare_product_distribution:
            comparison_results["comparison"]["product_distribution"] = _compare_product_distributions(
                comparison_results["individual_results"]
            )

        if request.compare_soa_yield:
            comparison_results["comparison"]["soa_yield"] = _compare_soa_yields(
                request.voc_names
            )

        if request.compare_radical_budget:
            comparison_results["comparison"]["radical_budget"] = _compare_radical_budgets(
                comparison_results["individual_results"]
            )

        # Generate comparison plots
        _generate_comparison_plots(comparison_results, output_dir)

        # Save comparison results
        with open(output_dir / "comparison_results.json", 'w') as f:
            json.dump(comparison_results, f, indent=2)

        job_manager.log_job(job_id, "Comparison workflow completed successfully")
        job_manager.update_job(
            job_id,
            status='completed',
            result_path=f"/data/output/{job_id}"
        )

    except Exception as e:
        job_manager.log_job(job_id, f"Comparison workflow failed: {str(e)}")
        job_manager.update_job(job_id, status='failed', error=str(e))


def _collect_mechanism_stats(output_dir: Path, voc_name: str) -> Dict[str, Any]:
    """Collect statistics about a mechanism."""
    stats = {
        "voc_name": voc_name,
        "num_species": 0,
        "num_reactions": 0,
        "num_radicals": 0,
        "num_stable": 0,
        "carbon_numbers": [],
        "functional_groups": {},
        "product_classes": {}
    }

    dictionary_file = output_dir / "dictionary.out"
    if dictionary_file.exists():
        with open(dictionary_file) as f:
            content = f.read()
            lines = [l.strip() for l in content.split('\n') if l.strip() and not l.startswith('!')]
            stats["num_species"] = len(lines)

            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    formula = parts[0]
                    # Count carbons
                    import re
                    c_match = re.search(r'C(\d+)', formula)
                    if c_match:
                        stats["carbon_numbers"].append(int(c_match.group(1)))

                    # Detect radicals (end in . or contain radical indicator)
                    if '.' in formula or formula.endswith('O') and 'O2' not in formula:
                        stats["num_radicals"] += 1
                    else:
                        stats["num_stable"] += 1

    mechanism_files = list(output_dir.glob("*.mec")) + list(output_dir.glob("*.k"))
    for mf in mechanism_files:
        with open(mf) as f:
            content = f.read()
            # Count reactions (lines with ->)
            stats["num_reactions"] += content.count('->')

    return stats


def _compare_mechanism_sizes(individual_results: Dict[str, Dict]) -> Dict[str, Any]:
    """Compare mechanism sizes across VOCs."""
    comparison = {
        "species_count": {},
        "reaction_count": {},
        "ranking_by_species": [],
        "ranking_by_reactions": []
    }

    for voc, stats in individual_results.items():
        if "error" not in stats:
            comparison["species_count"][voc] = stats.get("num_species", 0)
            comparison["reaction_count"][voc] = stats.get("num_reactions", 0)

    # Sort by size
    comparison["ranking_by_species"] = sorted(
        comparison["species_count"].items(), key=lambda x: x[1], reverse=True
    )
    comparison["ranking_by_reactions"] = sorted(
        comparison["reaction_count"].items(), key=lambda x: x[1], reverse=True
    )

    return comparison


def _compare_product_distributions(individual_results: Dict[str, Dict]) -> Dict[str, Any]:
    """Compare product distributions across VOCs."""
    comparison = {
        "carbon_distribution": {},
        "radical_fraction": {}
    }

    for voc, stats in individual_results.items():
        if "error" not in stats:
            carbon_nums = stats.get("carbon_numbers", [])
            if carbon_nums:
                comparison["carbon_distribution"][voc] = {
                    "min": min(carbon_nums),
                    "max": max(carbon_nums),
                    "mean": sum(carbon_nums) / len(carbon_nums)
                }

            total = stats.get("num_species", 0)
            radicals = stats.get("num_radicals", 0)
            if total > 0:
                comparison["radical_fraction"][voc] = radicals / total

    return comparison


def _compare_soa_yields(voc_names: List[str]) -> Dict[str, Any]:
    """
    Compare SOA yields from database.

    Note: SOA yields are typically determined experimentally and depend on
    NOx conditions, seed aerosol, etc. This function returns placeholder
    data - actual SOA yield comparison should be done via Box Model simulations.
    """
    comparison = {
        "high_nox_yields": {},
        "low_nox_yields": {},
        "ranking_high_nox": [],
        "ranking_low_nox": [],
        "note": "SOA yields require Box Model simulation for accurate comparison"
    }

    # SOA yield data is not stored in compound database
    # These would need to be calculated from Box Model runs
    for voc in voc_names:
        compound = compound_database.get_compound(voc)
        if compound:
            # Use category-based estimates as placeholders
            category = compound.category.lower() if compound.category else ""
            if "terpene" in category or "monoterpene" in category:
                comparison["high_nox_yields"][voc] = 0.15  # Typical monoterpene
                comparison["low_nox_yields"][voc] = 0.35
            elif "aromatic" in category:
                comparison["high_nox_yields"][voc] = 0.10  # Typical aromatic
                comparison["low_nox_yields"][voc] = 0.30
            elif "isoprene" in voc.lower():
                comparison["high_nox_yields"][voc] = 0.03  # Isoprene (low yield)
                comparison["low_nox_yields"][voc] = 0.10
            elif "alkane" in category:
                comparison["high_nox_yields"][voc] = 0.05  # Typical alkane
                comparison["low_nox_yields"][voc] = 0.15

    comparison["ranking_high_nox"] = sorted(
        comparison["high_nox_yields"].items(), key=lambda x: x[1], reverse=True
    )
    comparison["ranking_low_nox"] = sorted(
        comparison["low_nox_yields"].items(), key=lambda x: x[1], reverse=True
    )

    return comparison


def _compare_radical_budgets(individual_results: Dict[str, Dict]) -> Dict[str, Any]:
    """Compare radical budgets across VOCs."""
    comparison = {
        "radical_counts": {},
        "radical_fractions": {}
    }

    for voc, stats in individual_results.items():
        if "error" not in stats:
            comparison["radical_counts"][voc] = stats.get("num_radicals", 0)
            total = stats.get("num_species", 0)
            radicals = stats.get("num_radicals", 0)
            if total > 0:
                comparison["radical_fractions"][voc] = radicals / total

    return comparison


def _generate_comparison_plots(comparison_results: Dict, output_dir: Path):
    """Generate comparison visualization plots."""
    try:
        vocs = comparison_results["vocs"]
        individual = comparison_results["individual_results"]

        # Mechanism size comparison bar chart
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Species count
        ax = axes[0, 0]
        species_counts = [individual.get(v, {}).get("num_species", 0) for v in vocs]
        bars = ax.bar(range(len(vocs)), species_counts, color='steelblue')
        ax.set_xticks(range(len(vocs)))
        ax.set_xticklabels(vocs, rotation=45, ha='right')
        ax.set_ylabel('Number of Species')
        ax.set_title('Mechanism Size (Species)')
        for bar, count in zip(bars, species_counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                   str(count), ha='center', va='bottom', fontsize=8)

        # Reaction count
        ax = axes[0, 1]
        rxn_counts = [individual.get(v, {}).get("num_reactions", 0) for v in vocs]
        bars = ax.bar(range(len(vocs)), rxn_counts, color='coral')
        ax.set_xticks(range(len(vocs)))
        ax.set_xticklabels(vocs, rotation=45, ha='right')
        ax.set_ylabel('Number of Reactions')
        ax.set_title('Mechanism Size (Reactions)')
        for bar, count in zip(bars, rxn_counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                   str(count), ha='center', va='bottom', fontsize=8)

        # SOA yields (if available)
        ax = axes[1, 0]
        soa_data = comparison_results.get("comparison", {}).get("soa_yield", {})
        high_nox = [soa_data.get("high_nox_yields", {}).get(v, 0) for v in vocs]
        low_nox = [soa_data.get("low_nox_yields", {}).get(v, 0) for v in vocs]

        x = np.arange(len(vocs))
        width = 0.35
        ax.bar(x - width/2, high_nox, width, label='High NOx', color='orange')
        ax.bar(x + width/2, low_nox, width, label='Low NOx', color='green')
        ax.set_xticks(x)
        ax.set_xticklabels(vocs, rotation=45, ha='right')
        ax.set_ylabel('SOA Yield')
        ax.set_title('SOA Yields')
        ax.legend()

        # Radical fraction
        ax = axes[1, 1]
        radical_frac = [
            individual.get(v, {}).get("num_radicals", 0) /
            max(individual.get(v, {}).get("num_species", 1), 1)
            for v in vocs
        ]
        bars = ax.bar(range(len(vocs)), radical_frac, color='purple')
        ax.set_xticks(range(len(vocs)))
        ax.set_xticklabels(vocs, rotation=45, ha='right')
        ax.set_ylabel('Radical Fraction')
        ax.set_title('Radical Species Fraction')

        plt.tight_layout()
        plt.savefig(output_dir / "voc_comparison.png", dpi=150, bbox_inches='tight')
        plt.close()

    except Exception as e:
        logger.warning(f"Failed to generate comparison plots: {e}")


# ==============================================================================
# Extended API Endpoints - Mechanism Reduction (Recommendation #15)
# ==============================================================================

@app.post("/api/mechanism/reduce")
async def create_reduction_workflow(
    request: MechanismReductionRequest,
    background_tasks: BackgroundTasks
):
    """Create a mechanism reduction workflow."""
    source_job = job_manager.get_job(request.job_id)
    if not source_job:
        raise HTTPException(status_code=404, detail="Source job not found")
    if source_job.status != 'completed':
        raise HTTPException(status_code=400, detail="Source job not completed")

    job = job_manager.create_job("reduction", source_job.voc_name, {
        "source_job_id": request.job_id,
        "reduction_method": request.reduction_method,
        "target_species": request.target_species,
        "error_threshold": request.error_threshold,
        "min_species": request.min_species,
        "preserve_radicals": request.preserve_radicals,
        "preserve_soa_precursors": request.preserve_soa_precursors
    })

    background_tasks.add_task(run_reduction_workflow, job.id, request)
    return {"job_id": job.id, "status": "queued", "workflow": "reduction"}


async def run_reduction_workflow(job_id: str, request: MechanismReductionRequest):
    """Execute mechanism reduction workflow."""
    job_manager.log_job(job_id, f"Starting mechanism reduction using {request.reduction_method}")
    job_manager.update_job(job_id, status='running')

    try:
        source_dir = DATA_DIR / "output" / request.job_id
        output_dir = DATA_DIR / "output" / job_id
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load source mechanism - check multiple possible locations
        # Combined workflows store mechanism in a 'mechanism' subdirectory
        dictionary_file = None
        mechanism_files = []

        # Check direct location first (for generator-only jobs)
        if (source_dir / "dictionary.out").exists():
            dictionary_file = source_dir / "dictionary.out"
            mechanism_files = list(source_dir.glob("*.mec")) + list(source_dir.glob("*.k"))
        # Check combined workflow structure (mechanism subdirectory)
        elif (source_dir / "mechanism" / "dictionary.out").exists():
            dictionary_file = source_dir / "mechanism" / "dictionary.out"
            mech_subdir = source_dir / "mechanism"
            mechanism_files = list(mech_subdir.glob("*.mec")) + list(mech_subdir.glob("*.k"))
            job_manager.log_job(job_id, "Found mechanism in combined workflow subdirectory")
        # Check exports directory (another possible location)
        elif (source_dir / "exports" / "dictionary.out").exists():
            dictionary_file = source_dir / "exports" / "dictionary.out"
            exports_subdir = source_dir / "exports"
            mechanism_files = list(exports_subdir.glob("*.mec")) + list(exports_subdir.glob("*.k"))
            job_manager.log_job(job_id, "Found mechanism in exports subdirectory")

        if dictionary_file is None or not dictionary_file.exists():
            # List what files are available for debugging
            available_files = list(source_dir.rglob("*"))
            file_list = [str(f.relative_to(source_dir)) for f in available_files if f.is_file()][:20]
            raise FileNotFoundError(
                f"Source dictionary not found in job {request.job_id}. "
                f"Available files: {file_list}"
            )

        # Parse mechanism
        species = _parse_species_from_dictionary(dictionary_file)
        reactions = []
        for mf in mechanism_files:
            reactions.extend(_parse_reactions_from_mechanism(mf))

        # Validate that we have species and reactions to work with
        if not species:
            raise ValueError(
                f"No species found in dictionary file: {dictionary_file}. "
                "The file may be empty or in an unexpected format."
            )

        if not reactions:
            job_manager.log_job(job_id, "WARNING: No reactions found in mechanism files")
            # Continue anyway - some mechanisms may only have species definitions

        original_stats = {
            "num_species": len(species),
            "num_reactions": len(reactions)
        }
        job_manager.log_job(job_id, f"Original mechanism: {len(species)} species, {len(reactions)} reactions")

        # Apply reduction method
        if request.reduction_method == "drgep":
            reduced_species, reduced_reactions = _reduce_drgep(
                species, reactions, request.target_species,
                request.error_threshold, request.min_species,
                request.preserve_radicals, request.preserve_soa_precursors
            )
        elif request.reduction_method == "lumping":
            reduced_species, reduced_reactions = _reduce_lumping(
                species, reactions, request.min_species,
                request.preserve_radicals, request.preserve_soa_precursors
            )
        elif request.reduction_method == "pfa":
            reduced_species, reduced_reactions = _reduce_pfa(
                species, reactions, request.target_species,
                request.error_threshold, request.min_species
            )
        elif request.reduction_method == "sensitivity":
            reduced_species, reduced_reactions = _reduce_sensitivity(
                species, reactions, request.target_species,
                request.error_threshold, request.min_species
            )
        else:
            raise ValueError(f"Unknown reduction method: {request.reduction_method}")

        reduced_stats = {
            "num_species": len(reduced_species),
            "num_reactions": len(reduced_reactions),
            "reduction_ratio_species": (1 - len(reduced_species) / len(species)) if species else 0,
            "reduction_ratio_reactions": (1 - len(reduced_reactions) / len(reactions)) if reactions else 0
        }

        species_pct = f"{reduced_stats['reduction_ratio_species']:.1%}" if species else "N/A"
        reactions_pct = f"{reduced_stats['reduction_ratio_reactions']:.1%}" if reactions else "N/A"
        job_manager.log_job(job_id,
            f"Reduced mechanism: {len(reduced_species)} species ({species_pct} reduction), "
            f"{len(reduced_reactions)} reactions ({reactions_pct} reduction)"
        )

        # Write reduced mechanism files
        _write_reduced_dictionary(reduced_species, output_dir / "dictionary_reduced.out")
        _write_reduced_mechanism(reduced_reactions, output_dir / "mechanism_reduced.mec")

        # Save reduction report
        report = {
            "source_job_id": request.job_id,
            "reduction_method": request.reduction_method,
            "original_mechanism": original_stats,
            "reduced_mechanism": reduced_stats,
            "target_species": request.target_species,
            "parameters": {
                "error_threshold": request.error_threshold,
                "min_species": request.min_species,
                "preserve_radicals": request.preserve_radicals,
                "preserve_soa_precursors": request.preserve_soa_precursors
            },
            "removed_species": [s for s in species if s not in reduced_species]
        }

        with open(output_dir / "reduction_report.json", 'w') as f:
            json.dump(report, f, indent=2)

        job_manager.log_job(job_id, "Mechanism reduction completed successfully")
        job_manager.update_job(
            job_id,
            status='completed',
            result_path=f"/data/output/{job_id}"
        )

    except Exception as e:
        job_manager.log_job(job_id, f"Mechanism reduction failed: {str(e)}")
        job_manager.update_job(job_id, status='failed', error=str(e))


def _parse_species_from_dictionary(dictionary_file: Path) -> List[str]:
    """Parse species names from dictionary file."""
    species = []
    with open(dictionary_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('!'):
                parts = line.split()
                if parts:
                    species.append(parts[0])
    return species


def _parse_reactions_from_mechanism(mechanism_file: Path) -> List[Dict[str, Any]]:
    """Parse reactions from mechanism file."""
    reactions = []
    with open(mechanism_file) as f:
        for line in f:
            line = line.strip()
            if '->' in line:
                parts = line.split('->')
                if len(parts) == 2:
                    reactants = [r.strip() for r in parts[0].split('+')]
                    products_rate = parts[1].split(':')
                    products = [p.strip() for p in products_rate[0].split('+')]
                    rate = products_rate[1].strip() if len(products_rate) > 1 else ""
                    reactions.append({
                        "reactants": reactants,
                        "products": products,
                        "rate": rate,
                        "original": line
                    })
    return reactions


def _reduce_drgep(
    species: List[str],
    reactions: List[Dict],
    target_species: List[str],
    error_threshold: float,
    min_species: int,
    preserve_radicals: bool,
    preserve_soa_precursors: bool
) -> tuple:
    """
    Directed Relation Graph with Error Propagation (DRGEP) reduction.

    Based on Pepiot-Desjardins & Pitsch (2008) "An efficient error-propagation-based
    reduction method for large chemical kinetic mechanisms"
    """
    # Build species interaction graph
    graph = {s: {} for s in species}

    for rxn in reactions:
        all_species = set(rxn["reactants"] + rxn["products"])
        for s1 in all_species:
            if s1 in graph:
                for s2 in all_species:
                    if s2 != s1 and s2 in graph:
                        # Simple interaction coefficient (can be refined with rate analysis)
                        graph[s1][s2] = graph[s1].get(s2, 0) + 1

    # Normalize interaction coefficients
    for s1, connections in graph.items():
        total = sum(connections.values()) if connections else 1
        for s2 in connections:
            graph[s1][s2] /= total

    # Auto-detect target species if not provided
    if not target_species:
        # Use parent VOC and key oxidation products
        target_species = [s for s in species if len(s) >= 4 and not s.endswith('.')][:5]

    # Essential species to preserve
    essential = set(target_species)

    if preserve_radicals:
        essential.update(s for s in species if s.endswith('.') or 'O2' in s)

    if preserve_soa_precursors:
        # Preserve low-volatility species (typically have many O atoms)
        essential.update(s for s in species if s.count('O') >= 4 and not s.endswith('.'))

    # DRGEP path-based analysis
    def drgep_coefficient(source: str, target: str, visited: set = None) -> float:
        if visited is None:
            visited = set()
        if source == target:
            return 1.0
        if source in visited:
            return 0.0

        visited.add(source)
        max_coeff = 0.0

        for neighbor, weight in graph.get(source, {}).items():
            if neighbor not in visited:
                path_coeff = weight * drgep_coefficient(neighbor, target, visited.copy())
                max_coeff = max(max_coeff, path_coeff)

        return max_coeff

    # Calculate importance of each species
    importance = {}
    for s in species:
        max_importance = 0.0
        for target in target_species:
            if target in species:
                coeff = drgep_coefficient(s, target)
                max_importance = max(max_importance, coeff)
        importance[s] = max_importance

    # Keep species above threshold or in essential set
    threshold = error_threshold
    kept_species = set()
    for s, imp in sorted(importance.items(), key=lambda x: x[1], reverse=True):
        if s in essential or imp >= threshold:
            kept_species.add(s)
        if len(kept_species) >= min_species and len(kept_species) >= len(essential):
            break

    # Ensure minimum species
    while len(kept_species) < min_species and len(kept_species) < len(species):
        remaining = [(s, imp) for s, imp in importance.items() if s not in kept_species]
        if remaining:
            remaining.sort(key=lambda x: x[1], reverse=True)
            kept_species.add(remaining[0][0])
        else:
            break

    # Filter reactions to only include those with kept species
    kept_reactions = []
    for rxn in reactions:
        all_species = set(rxn["reactants"] + rxn["products"])
        if all_species.issubset(kept_species):
            kept_reactions.append(rxn)

    return list(kept_species), kept_reactions


def _reduce_lumping(
    species: List[str],
    reactions: List[Dict],
    min_species: int,
    preserve_radicals: bool,
    preserve_soa_precursors: bool
) -> tuple:
    """
    Species lumping reduction.

    Groups similar species based on structural similarity and replaces
    with representative species.
    """
    # Group species by carbon number and functional groups
    groups = {}

    for s in species:
        # Extract carbon number
        import re
        c_match = re.search(r'C(\d+)', s)
        c_num = int(c_match.group(1)) if c_match else 0

        # Count functional groups
        o_count = s.count('O')
        n_count = s.count('N')
        is_radical = s.endswith('.') or ('O2' in s and not s.startswith('GO'))

        key = (c_num, o_count // 2, n_count, is_radical)  # Bin O atoms
        if key not in groups:
            groups[key] = []
        groups[key].append(s)

    # Select representative from each group
    kept_species = set()
    lumping_map = {}  # Original -> Representative

    for key, group in groups.items():
        c_num, o_bin, n_count, is_radical = key

        # Keep radicals if requested
        if is_radical and preserve_radicals:
            for s in group:
                kept_species.add(s)
                lumping_map[s] = s
            continue

        # Keep low-volatility species if requested
        if preserve_soa_precursors and o_bin >= 2:
            for s in group:
                kept_species.add(s)
                lumping_map[s] = s
            continue

        # Otherwise, select representative (first or most common)
        representative = group[0]
        kept_species.add(representative)
        for s in group:
            lumping_map[s] = representative

    # Ensure minimum species
    if len(kept_species) < min_species:
        remaining = [s for s in species if s not in kept_species]
        kept_species.update(remaining[:min_species - len(kept_species)])
        for s in kept_species:
            lumping_map[s] = s

    # Update reactions with lumped species
    kept_reactions = []
    seen_reactions = set()

    for rxn in reactions:
        new_reactants = [lumping_map.get(r, r) for r in rxn["reactants"]]
        new_products = [lumping_map.get(p, p) for p in rxn["products"]]

        # Check if all species are in kept set
        all_species = set(new_reactants + new_products)
        if not all_species.issubset(kept_species):
            continue

        # Create unique reaction key to avoid duplicates
        rxn_key = (tuple(sorted(new_reactants)), tuple(sorted(new_products)))
        if rxn_key not in seen_reactions:
            seen_reactions.add(rxn_key)
            kept_reactions.append({
                "reactants": new_reactants,
                "products": new_products,
                "rate": rxn["rate"],
                "original": f"{' + '.join(new_reactants)} -> {' + '.join(new_products)}"
            })

    return list(kept_species), kept_reactions


def _reduce_pfa(
    species: List[str],
    reactions: List[Dict],
    target_species: List[str],
    error_threshold: float,
    min_species: int
) -> tuple:
    """
    Path Flux Analysis (PFA) reduction.

    Based on Sun et al. (2010) "Path flux analysis for the reduction of
    detailed chemical kinetic mechanisms"
    """
    # Build flux graph
    production = {s: [] for s in species}
    consumption = {s: [] for s in species}

    for i, rxn in enumerate(reactions):
        for r in rxn["reactants"]:
            if r in consumption:
                consumption[r].append(i)
        for p in rxn["products"]:
            if p in production:
                production[p].append(i)

    # Calculate path flux coefficients
    def path_flux(source: str, target: str, max_depth: int = 5) -> float:
        if source == target:
            return 1.0
        if max_depth <= 0:
            return 0.0

        total_flux = 0.0
        for rxn_idx in production.get(source, []):
            rxn = reactions[rxn_idx]
            for next_species in rxn["reactants"]:
                if next_species != source:
                    total_flux += path_flux(next_species, target, max_depth - 1)

        return min(total_flux, 1.0)

    # Auto-detect targets if not provided
    if not target_species:
        target_species = [s for s in species if len(s) >= 4][:5]

    # Calculate importance
    importance = {}
    for s in species:
        max_flux = 0.0
        for target in target_species:
            if target in species:
                flux = path_flux(s, target)
                max_flux = max(max_flux, flux)
        importance[s] = max_flux

    # Keep species above threshold
    kept_species = {s for s, imp in importance.items() if imp >= error_threshold}
    kept_species.update(target_species)

    while len(kept_species) < min_species and len(kept_species) < len(species):
        remaining = [(s, imp) for s, imp in importance.items() if s not in kept_species]
        if remaining:
            remaining.sort(key=lambda x: x[1], reverse=True)
            kept_species.add(remaining[0][0])
        else:
            break

    # Filter reactions
    kept_reactions = [
        rxn for rxn in reactions
        if set(rxn["reactants"] + rxn["products"]).issubset(kept_species)
    ]

    return list(kept_species), kept_reactions


def _reduce_sensitivity(
    species: List[str],
    reactions: List[Dict],
    target_species: List[str],
    error_threshold: float,
    min_species: int
) -> tuple:
    """
    Sensitivity-based reduction.

    Removes species/reactions with low sensitivity coefficients.
    """
    # Simple connectivity-based sensitivity
    connectivity = {s: 0 for s in species}

    for rxn in reactions:
        for s in rxn["reactants"] + rxn["products"]:
            if s in connectivity:
                connectivity[s] += 1

    # Normalize
    max_conn = max(connectivity.values()) if connectivity else 1
    sensitivity = {s: c / max_conn for s, c in connectivity.items()}

    # Auto-detect targets
    if not target_species:
        target_species = [s for s in species if len(s) >= 4][:5]

    # Keep high-sensitivity species and targets
    kept_species = {s for s, sens in sensitivity.items() if sens >= error_threshold}
    kept_species.update(target_species)

    while len(kept_species) < min_species and len(kept_species) < len(species):
        remaining = [(s, sens) for s, sens in sensitivity.items() if s not in kept_species]
        if remaining:
            remaining.sort(key=lambda x: x[1], reverse=True)
            kept_species.add(remaining[0][0])
        else:
            break

    kept_reactions = [
        rxn for rxn in reactions
        if set(rxn["reactants"] + rxn["products"]).issubset(kept_species)
    ]

    return list(kept_species), kept_reactions


def _write_reduced_dictionary(species: List[str], output_path: Path):
    """Write reduced species dictionary."""
    with open(output_path, 'w') as f:
        f.write("! Reduced mechanism dictionary\n")
        f.write(f"! {len(species)} species\n")
        f.write("!\n")
        for s in sorted(species):
            f.write(f"{s}\n")


def _write_reduced_mechanism(reactions: List[Dict], output_path: Path):
    """Write reduced mechanism file."""
    with open(output_path, 'w') as f:
        f.write("! Reduced mechanism\n")
        f.write(f"! {len(reactions)} reactions\n")
        f.write("!\n")
        for rxn in reactions:
            reactants = " + ".join(rxn["reactants"])
            products = " + ".join(rxn["products"])
            rate = rxn.get("rate", "")
            if rate:
                f.write(f"{reactants} -> {products} : {rate}\n")
            else:
                f.write(f"{reactants} -> {products}\n")


# ==============================================================================
# Extended API Endpoints - 3D Structure Visualization (Recommendation #13)
# ==============================================================================

def _sanitize_smiles_for_3d(smiles: str) -> str:
    """
    Sanitize SMILES for 3D embedding by handling radicals and other problematic features.

    RDKit's embedding can fail for:
    - Radical species (e.g., [O], [N], etc.)
    - Some charged species
    - Unusual valence states

    This function attempts to convert these to embeddable forms while preserving
    the overall structure for visualization.
    """
    import re

    # Common radical/problematic patterns and their substitutions for 3D visualization
    # These substitutions maintain the basic structure while making embedding possible
    # Order matters - more specific patterns first
    substitutions = [
        # Peroxynitrate with terminal radical: [N+](=O)[O-]O[O] -> [N+](=O)[O-]O
        (r'\[N\+\]\(=O\)\[O\-\]O\[O\]', '[N+](=O)[O-]O'),
        # Terminal peroxy radical O[O] -> OO (peroxide for visualization)
        (r'O\[O\]$', 'OO'),
        # Internal peroxy radical O[O] -> OO
        (r'O\[O\]', 'OO'),
        # Terminal oxygen radical [O] -> O
        (r'\[O\]$', 'O'),
        # Oxygen radical not adjacent to negative charge: [O] -> O
        (r'\[O\](?!\-)', 'O'),
        # Nitrogen radicals (not charged)
        (r'\[N\](?!\+|\-)', 'N'),
        # Carbon radicals [CH2] -> C (methylene)
        (r'\[CH2\]', 'C'),
        # Carbon radicals [CH] -> C
        (r'\[CH\]', 'C'),
        # Remove any remaining bare radical markers
        (r'\[\.\]', ''),
    ]

    sanitized = smiles
    for pattern, replacement in substitutions:
        sanitized = re.sub(pattern, replacement, sanitized)

    return sanitized


@app.get("/api/structure/{compound_name}/3d")
async def get_3d_structure(compound_name: str, format: str = "mol"):
    """
    Get 3D molecular structure for interactive visualization.

    Returns MOL/SDF format with 3D coordinates for use with 3Dmol.js or similar.

    Handles difficult structures including:
    - Radical species (oxygen, nitrogen radicals)
    - Peroxy compounds
    - Complex nitrate/peroxy combinations
    """
    # Get SMILES for compound
    smiles = None

    # Check if compound_name is already a SMILES string
    if any(c in compound_name for c in ['=', '(', ')', '[', ']', '#', '@']):
        # Looks like SMILES notation
        smiles = compound_name
    else:
        # Try compound database first
        compound = compound_database.get_compound(compound_name)
        if compound:
            smiles = compound.smiles
        else:
            # Try reaction tree conversion (GECKO code -> SMILES)
            gecko_formula = compound_name
            converted = reaction_tree.gecko_to_smiles(gecko_formula)
            if converted != gecko_formula:
                smiles = converted

    if not smiles:
        raise HTTPException(
            status_code=404,
            detail=f"Cannot find SMILES for '{compound_name}'. Try using the compound's SMILES directly."
        )

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # First try with original SMILES
        original_smiles = smiles
        mol = Chem.MolFromSmiles(smiles)

        # If original fails to parse, try sanitizing
        if mol is None:
            smiles = _sanitize_smiles_for_3d(smiles)
            mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid SMILES: {original_smiles}. "
                       "This structure cannot be parsed by RDKit."
            )

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)

        # Try multiple embedding strategies
        result = -1

        # Strategy 1: Standard embedding
        if result == -1:
            try:
                result = AllChem.EmbedMolecule(mol, randomSeed=42)
            except Exception:
                pass

        # Strategy 2: Random coordinates
        if result == -1:
            try:
                result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
            except Exception:
                pass

        # Strategy 3: ETKDG with more permissive settings
        if result == -1:
            try:
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.useRandomCoords = True
                params.maxIterations = 500
                result = AllChem.EmbedMolecule(mol, params)
            except Exception:
                pass

        # Strategy 4: Try with sanitized SMILES if original didn't work
        if result == -1 and smiles == original_smiles:
            sanitized = _sanitize_smiles_for_3d(smiles)
            if sanitized != smiles:
                try:
                    mol = Chem.MolFromSmiles(sanitized)
                    if mol:
                        mol = Chem.AddHs(mol)
                        result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
                        if result == 0:
                            smiles = sanitized  # Use sanitized for response
                except Exception:
                    pass

        # Strategy 5: Generate 2D coordinates and project to 3D
        if result == -1:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mol = Chem.AddHs(mol)
                    AllChem.Compute2DCoords(mol)
                    # Create fake 3D by adding z=0 to all atoms
                    conf = mol.GetConformer()
                    for i in range(mol.GetNumAtoms()):
                        pos = conf.GetAtomPosition(i)
                        conf.SetAtomPosition(i, (pos.x, pos.y, 0.0))
                    result = 0  # Mark as success (2D projected)
            except Exception:
                pass

        if result == -1:
            # All embedding attempts failed
            raise HTTPException(
                status_code=400,
                detail=f"No 3D structure for: {original_smiles}. "
                       "This may be a radical species or unusual structure that cannot be embedded. "
                       "Try viewing a related stable compound instead."
            )

        # Optimize geometry if embedding succeeded
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            # MMFF optimization can fail for some structures, try UFF as fallback
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                # Continue without optimization if both fail
                pass

        if format == "mol":
            mol_block = Chem.MolToMolBlock(mol)
            return {"format": "mol", "data": mol_block, "smiles": smiles}
        elif format == "sdf":
            sdf_block = Chem.MolToMolBlock(mol)
            return {"format": "sdf", "data": sdf_block, "smiles": smiles}
        elif format == "xyz":
            # Generate XYZ format
            conf = mol.GetConformer()
            xyz_lines = [str(mol.GetNumAtoms()), compound_name]
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                xyz_lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
            return {"format": "xyz", "data": "\n".join(xyz_lines), "smiles": smiles}
        elif format == "pdb":
            pdb_block = Chem.MolToPDBBlock(mol)
            return {"format": "pdb", "data": pdb_block, "smiles": smiles}
        else:
            raise HTTPException(status_code=400, detail=f"Unknown format: {format}")

    except ImportError:
        raise HTTPException(status_code=501, detail="RDKit not available for 3D structure generation")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"3D structure generation failed: {str(e)}")


@app.get("/api/jobs/{job_id}/structures/3d")
async def get_job_3d_structures(job_id: str, max_structures: int = 20):
    """Get 3D structures for all species in a completed job."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != 'completed':
        raise HTTPException(status_code=400, detail="Job not completed")

    output_dir = DATA_DIR / "output" / job_id
    tree_file = output_dir / "reaction_tree.json"

    if not tree_file.exists():
        raise HTTPException(status_code=404, detail="Reaction tree not found")

    with open(tree_file) as f:
        tree_data = json.load(f)

    structures = []
    nodes = tree_data.get("nodes", [])[:max_structures]

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        for node in nodes:
            smiles = node.get("smiles")
            if not smiles or smiles == node.get("id"):
                continue

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue

                mol = Chem.AddHs(mol)
                result = AllChem.EmbedMolecule(mol, randomSeed=42)
                if result == -1:
                    result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
                if result == 0:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=100)

                mol_block = Chem.MolToMolBlock(mol)
                structures.append({
                    "id": node.get("id"),
                    "smiles": smiles,
                    "label": node.get("label", ""),
                    "mol_block": mol_block
                })
            except Exception:
                continue

        return {"job_id": job_id, "structures": structures, "count": len(structures)}

    except ImportError:
        raise HTTPException(status_code=501, detail="RDKit not available for 3D structure generation")


# ==============================================================================
# Extended API Endpoints - Enhanced Visualization
# ==============================================================================

@app.post("/api/jobs/{job_id}/visualize")
async def generate_enhanced_visualization(
    job_id: str,
    layout: str = "hierarchical",
    show_branching: bool = True,
    color_by_reaction: bool = True,
    format: str = "png"
):
    """Generate enhanced pathway visualization with custom options."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != 'completed':
        raise HTTPException(status_code=400, detail="Job not completed")

    output_dir = DATA_DIR / "output" / job_id
    tree_file = output_dir / "reaction_tree.json"

    if not tree_file.exists():
        raise HTTPException(status_code=404, detail="Reaction tree not found")

    with open(tree_file) as f:
        tree_data = json.load(f)

    # Check if graphviz is available
    if not enhanced_visualizer.HAS_GRAPHVIZ:
        raise HTTPException(
            status_code=501,
            detail="Diagram regeneration requires graphviz. Install with: brew install graphviz (Mac) or apt-get install graphviz (Linux)"
        )

    try:
        # Create visualizer without passing tree_data to __init__
        # tree_data must be loaded via from_tree_data() method
        visualizer = enhanced_visualizer.EnhancedPathwayVisualizer()
        visualizer.from_tree_data(tree_data)

        output_name = f"pathway_{layout}"
        output_path = output_dir / output_name

        visualizer.generate_diagram(
            output_path=str(output_path),
            title=f"Oxidation Mechanism: {job.voc_name.upper()}",
            layout=layout,
            show_branching_ratios=show_branching,
            color_by_reaction_type=color_by_reaction,
            format=format
        )

        return {
            "success": True,
            "diagram_url": f"/data/output/{job_id}/{output_name}.{format}",
            "layout": layout,
            "format": format
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Visualization failed: {str(e)}")


# ==============================================================================
# Startup Events
# ==============================================================================

@app.on_event("startup")
async def startup_event():
    """Cleanup old jobs on startup."""
    job_manager.cleanup_old_jobs()
    logger.info("GECKO-A Web Interface started")
    logger.info(f"Compound database: {len(compound_database.COMPOUND_DATABASE)} compounds available")


# ==============================================================================
# Main Entry Point
# ==============================================================================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=5005, log_level="info")
