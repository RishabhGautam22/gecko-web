"""
GECKO-A Combined Workflow Module

This module provides a unified workflow that combines:
1. GECKO-A Mechanism Generator
2. BOXMODEL4GECKO Simulation
3. Mass Balance Verification
4. Post-processing and Analysis
5. Diagram Generation

The workflow runs all steps sequentially with proper error handling,
produces structured output directories, and generates comprehensive reports.

Author: GECKO-A Development Team
"""

import os
import json
import shutil
import logging
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field, asdict
from enum import Enum

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ==============================================================================
# Workflow Configuration
# ==============================================================================

class WorkflowStage(Enum):
    """Stages in the combined workflow."""
    INITIALIZATION = "initialization"
    MECHANISM_GENERATION = "mechanism_generation"
    MASS_BALANCE_CHECK = "mass_balance_check"
    BOX_MODEL_SIMULATION = "box_model_simulation"
    POST_PROCESSING = "post_processing"
    DIAGRAM_GENERATION = "diagram_generation"
    REPORT_GENERATION = "report_generation"
    FINALIZATION = "finalization"


class WorkflowStatus(Enum):
    """Status of workflow execution."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"
    WARNING = "warning"


@dataclass
class StageResult:
    """Result of a single workflow stage."""
    stage: WorkflowStage
    status: WorkflowStatus
    start_time: float
    end_time: float
    duration_seconds: float
    message: str = ""
    error: Optional[str] = None
    outputs: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            'stage': self.stage.value,
            'status': self.status.value,
            'start_time': datetime.fromtimestamp(self.start_time).isoformat(),
            'end_time': datetime.fromtimestamp(self.end_time).isoformat(),
            'duration_seconds': self.duration_seconds,
            'message': self.message,
            'error': self.error,
            'outputs': self.outputs,
            'warnings': self.warnings
        }


@dataclass
class WorkflowResult:
    """Complete result of workflow execution."""
    workflow_id: str
    voc_name: str
    status: WorkflowStatus
    start_time: float
    end_time: float
    total_duration_seconds: float
    stages: List[StageResult]
    output_directory: str
    summary: Dict[str, Any] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        """Return True if workflow completed successfully."""
        return self.status in (WorkflowStatus.COMPLETED, WorkflowStatus.WARNING)

    @property
    def error(self) -> Optional[str]:
        """Return error message from first failed stage, if any."""
        for stage in self.stages:
            if stage.status == WorkflowStatus.FAILED and stage.error:
                return f"Stage '{stage.stage.value}' failed: {stage.error}"
        if self.status == WorkflowStatus.FAILED:
            return "Workflow failed (unknown error)"
        return None

    def to_dict(self) -> Dict:
        return {
            'workflow_id': self.workflow_id,
            'voc_name': self.voc_name,
            'status': self.status.value,
            'start_time': datetime.fromtimestamp(self.start_time).isoformat(),
            'end_time': datetime.fromtimestamp(self.end_time).isoformat(),
            'total_duration_seconds': self.total_duration_seconds,
            'stages': [s.to_dict() for s in self.stages],
            'output_directory': self.output_directory,
            'summary': self.summary
        }


@dataclass
class WorkflowConfig:
    """Configuration for combined workflow execution."""
    # VOC identification
    voc_name: str
    voc_formula: Optional[str] = None

    # Output directory
    output_base: Optional[Path] = None

    # Workflow stage control
    run_generator: bool = True
    run_boxmodel: bool = True
    run_postprocessing: bool = True
    verify_mass_balance: bool = True

    # Generator parameters
    vapor_pressure_threshold: float = -13.0
    max_generations: int = 2
    output_format: str = "kpp"  # kpp, mcm, facsimile
    generator_options: Optional[Dict[str, Any]] = None

    # Box model parameters
    run_box_model: bool = True  # Alias for run_boxmodel
    temperature_k: float = 298.0
    rh_percent: float = 50.0
    latitude_degrees: float = 45.0
    simulation_hours: float = 24.0
    initial_o3_ppb: float = 40.0
    initial_nox_ppb: float = 10.0
    seed_aerosol_ug_m3: float = 10.0
    dilution_rate_s1: float = 0.0
    boxmodel_options: Optional[Dict[str, Any]] = None

    # Verification options
    run_mass_balance: bool = True
    mass_balance_tolerance: float = 0.05

    # Output options
    generate_diagrams: bool = True
    generate_pdf_report: bool = True
    generate_all_formats: bool = True
    max_diagram_nodes: int = 50
    diagram_max_depth: int = 3

    # Processing options
    skip_on_cache_hit: bool = True
    preserve_intermediate: bool = False
    verbose: bool = False


# ==============================================================================
# Output Directory Structure
# ==============================================================================

class OutputStructure:
    """
    Manages structured output directory layout.

    Directory structure:
    output_dir/
    ├── mechanism/
    │   ├── dictionary.out
    │   ├── reactions.mec
    │   ├── species.kpp
    │   └── ...
    ├── simulation/
    │   ├── concentrations.nc
    │   ├── time_series.csv
    │   └── ...
    ├── analysis/
    │   ├── aerosol_data.csv
    │   ├── mass_balance.json
    │   ├── partitioning_summary.json
    │   └── ...
    ├── plots/
    │   ├── volatility_distribution.png
    │   ├── soa_yield.png
    │   ├── van_krevelen.png
    │   └── ...
    ├── diagrams/
    │   ├── reaction_tree.json
    │   ├── pathway_diagram.png
    │   ├── pathway_diagram.svg
    │   └── mechanism_diagram.png
    ├── exports/
    │   ├── mechanism.kpp
    │   ├── mechanism.mcm
    │   └── mechanism.fac
    ├── reports/
    │   ├── workflow_report.json
    │   ├── summary.txt
    │   └── report.pdf (if enabled)
    └── logs/
        └── workflow.log
    """

    SUBDIRS = [
        'mechanism',
        'simulation',
        'analysis',
        'plots',
        'diagrams',
        'exports',
        'reports',
        'logs'
    ]

    def __init__(self, base_dir: Path, workflow_id: str = None):
        self.base_dir = Path(base_dir)
        self.workflow_id = workflow_id
        # Use base_dir directly as root (don't create nested workflow_id folder)
        # since main.py already creates job-specific directory
        self.root = self.base_dir
        self._create_structure()

    def _create_structure(self):
        """Create the directory structure."""
        self.root.mkdir(parents=True, exist_ok=True)
        for subdir in self.SUBDIRS:
            (self.root / subdir).mkdir(exist_ok=True)

    @property
    def mechanism_dir(self) -> Path:
        return self.root / "mechanism"

    @property
    def simulation_dir(self) -> Path:
        return self.root / "simulation"

    @property
    def analysis_dir(self) -> Path:
        return self.root / "analysis"

    @property
    def plots_dir(self) -> Path:
        return self.root / "plots"

    @property
    def diagrams_dir(self) -> Path:
        return self.root / "diagrams"

    @property
    def exports_dir(self) -> Path:
        return self.root / "exports"

    @property
    def reports_dir(self) -> Path:
        return self.root / "reports"

    @property
    def logs_dir(self) -> Path:
        return self.root / "logs"

    def get_log_path(self) -> Path:
        return self.logs_dir / "workflow.log"

    def list_outputs(self) -> Dict[str, List[str]]:
        """List all output files by category."""
        outputs = {}
        for subdir in self.SUBDIRS:
            subdir_path = self.root / subdir
            if subdir_path.exists():
                outputs[subdir] = [f.name for f in subdir_path.iterdir() if f.is_file()]
        return outputs


# ==============================================================================
# Combined Workflow Executor
# ==============================================================================

class CombinedWorkflow:
    """
    Executes the combined GECKO-A + Box Model workflow.

    This class orchestrates:
    1. Mechanism generation with GECKO-A
    2. Mass balance verification
    3. Box model simulation
    4. Post-processing and analysis
    5. Diagram and report generation

    Usage:
        config = WorkflowConfig(voc_name="alpha-pinene")
        workflow = CombinedWorkflow(config, output_dir="/path/to/output")
        result = workflow.run()
    """

    def __init__(self,
                 config: WorkflowConfig,
                 output_dir: Optional[Path] = None,
                 workflow_id: Optional[str] = None,
                 progress_callback: Optional[Callable[[str, float], None]] = None):
        """
        Initialize workflow.

        Args:
            config: Workflow configuration
            output_dir: Base output directory (can also be set via config.output_base)
            workflow_id: Optional workflow ID (auto-generated if not provided)
            progress_callback: Optional callback for progress updates (message, fraction)
        """
        self.config = config

        # Get output directory from config or parameter
        if output_dir:
            self.output_dir = Path(output_dir)
        elif config.output_base:
            self.output_dir = Path(config.output_base)
        else:
            raise ValueError("output_dir must be provided either as parameter or in config.output_base")

        self.workflow_id = workflow_id or self._generate_id()
        self.progress_callback = progress_callback

        # Initialize output structure
        self.output = OutputStructure(self.output_dir, self.workflow_id)

        # Setup logging
        self._setup_logging()

        # Track results
        self.stages: List[StageResult] = []
        self._start_time: float = 0
        self._current_stage: Optional[WorkflowStage] = None

    def _generate_id(self) -> str:
        """Generate a unique workflow ID."""
        import uuid
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        short_uuid = uuid.uuid4().hex[:8]
        voc_clean = self.config.voc_name.lower().replace('-', '').replace(' ', '')[:10]
        return f"wf_{voc_clean}_{timestamp}_{short_uuid}"

    def _setup_logging(self):
        """Configure logging for the workflow."""
        log_file = self.output.get_log_path()

        # File handler for detailed logs
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))

        # Add handler to logger
        self._logger = logging.getLogger(f"workflow.{self.workflow_id}")
        self._logger.addHandler(file_handler)
        self._logger.setLevel(logging.DEBUG)

    def _log(self, message: str, level: str = "info"):
        """Log a message and optionally update progress."""
        log_func = getattr(self._logger, level, self._logger.info)
        log_func(message)

        if self.config.verbose:
            print(f"[{self.workflow_id}] {message}")

    def _update_progress(self, message: str, fraction: float):
        """Update progress if callback is provided."""
        if self.progress_callback:
            self.progress_callback(message, fraction)

    def _run_stage(self, stage: WorkflowStage, func: Callable[[], Dict],
                   progress_start: float, progress_end: float) -> StageResult:
        """Execute a workflow stage with timing and error handling."""
        self._current_stage = stage
        start_time = time.time()

        self._log(f"Starting stage: {stage.value}")
        self._update_progress(f"Running {stage.value}...", progress_start)

        try:
            result_data = func()

            end_time = time.time()
            duration = end_time - start_time

            status = WorkflowStatus.COMPLETED
            if result_data.get('warnings'):
                status = WorkflowStatus.WARNING

            result = StageResult(
                stage=stage,
                status=status,
                start_time=start_time,
                end_time=end_time,
                duration_seconds=duration,
                message=result_data.get('message', f'{stage.value} completed'),
                outputs=result_data.get('outputs', {}),
                warnings=result_data.get('warnings', [])
            )

            self._log(f"Completed stage {stage.value} in {duration:.2f}s")
            self._update_progress(f"Completed {stage.value}", progress_end)

            return result

        except Exception as e:
            end_time = time.time()
            duration = end_time - start_time

            self._log(f"Stage {stage.value} failed: {e}", level="error")

            return StageResult(
                stage=stage,
                status=WorkflowStatus.FAILED,
                start_time=start_time,
                end_time=end_time,
                duration_seconds=duration,
                message=f"Stage failed: {stage.value}",
                error=str(e)
            )

    def run(self) -> WorkflowResult:
        """
        Execute the complete workflow.

        Returns:
            WorkflowResult with execution details and outputs
        """
        self._start_time = time.time()
        self._log(f"Starting combined workflow for {self.config.voc_name}")

        # Stage progress ranges (start, end)
        progress_map = {
            WorkflowStage.INITIALIZATION: (0.0, 0.05),
            WorkflowStage.MECHANISM_GENERATION: (0.05, 0.35),
            WorkflowStage.MASS_BALANCE_CHECK: (0.35, 0.45),
            WorkflowStage.BOX_MODEL_SIMULATION: (0.45, 0.65),
            WorkflowStage.POST_PROCESSING: (0.65, 0.75),
            WorkflowStage.DIAGRAM_GENERATION: (0.75, 0.85),
            WorkflowStage.REPORT_GENERATION: (0.85, 0.95),
            WorkflowStage.FINALIZATION: (0.95, 1.0),
        }

        workflow_failed = False

        # Stage 1: Initialization
        result = self._run_stage(
            WorkflowStage.INITIALIZATION,
            self._stage_initialization,
            *progress_map[WorkflowStage.INITIALIZATION]
        )
        self.stages.append(result)
        if result.status == WorkflowStatus.FAILED:
            workflow_failed = True

        # Stage 2: Mechanism Generation
        if not workflow_failed:
            result = self._run_stage(
                WorkflowStage.MECHANISM_GENERATION,
                self._stage_mechanism_generation,
                *progress_map[WorkflowStage.MECHANISM_GENERATION]
            )
            self.stages.append(result)
            if result.status == WorkflowStatus.FAILED:
                workflow_failed = True

        # Stage 3: Mass Balance Check
        if not workflow_failed and self.config.run_mass_balance:
            result = self._run_stage(
                WorkflowStage.MASS_BALANCE_CHECK,
                self._stage_mass_balance,
                *progress_map[WorkflowStage.MASS_BALANCE_CHECK]
            )
            self.stages.append(result)
            # Mass balance failure is a warning, not a workflow failure

        # Stage 4: Box Model Simulation
        if not workflow_failed and self.config.run_box_model:
            result = self._run_stage(
                WorkflowStage.BOX_MODEL_SIMULATION,
                self._stage_box_model,
                *progress_map[WorkflowStage.BOX_MODEL_SIMULATION]
            )
            self.stages.append(result)
            if result.status == WorkflowStatus.FAILED:
                workflow_failed = True

        # Stage 5: Post-processing
        if not workflow_failed:
            result = self._run_stage(
                WorkflowStage.POST_PROCESSING,
                self._stage_post_processing,
                *progress_map[WorkflowStage.POST_PROCESSING]
            )
            self.stages.append(result)
            # Post-processing failure is a warning

        # Stage 6: Diagram Generation
        if not workflow_failed and self.config.generate_diagrams:
            result = self._run_stage(
                WorkflowStage.DIAGRAM_GENERATION,
                self._stage_diagram_generation,
                *progress_map[WorkflowStage.DIAGRAM_GENERATION]
            )
            self.stages.append(result)
            # Diagram failure is a warning

        # Stage 7: Report Generation
        if not workflow_failed:
            result = self._run_stage(
                WorkflowStage.REPORT_GENERATION,
                self._stage_report_generation,
                *progress_map[WorkflowStage.REPORT_GENERATION]
            )
            self.stages.append(result)

        # Stage 8: Finalization
        result = self._run_stage(
            WorkflowStage.FINALIZATION,
            self._stage_finalization,
            *progress_map[WorkflowStage.FINALIZATION]
        )
        self.stages.append(result)

        # Build final result
        end_time = time.time()
        total_duration = end_time - self._start_time

        overall_status = WorkflowStatus.COMPLETED
        if workflow_failed:
            overall_status = WorkflowStatus.FAILED
        elif any(s.status == WorkflowStatus.WARNING for s in self.stages):
            overall_status = WorkflowStatus.WARNING

        workflow_result = WorkflowResult(
            workflow_id=self.workflow_id,
            voc_name=self.config.voc_name,
            status=overall_status,
            start_time=self._start_time,
            end_time=end_time,
            total_duration_seconds=total_duration,
            stages=self.stages,
            output_directory=str(self.output.root),
            summary=self._build_summary()
        )

        # Save workflow result
        result_file = self.output.reports_dir / "workflow_result.json"
        with open(result_file, 'w') as f:
            json.dump(workflow_result.to_dict(), f, indent=2)

        self._log(f"Workflow completed with status: {overall_status.value}")
        self._update_progress("Workflow completed", 1.0)

        return workflow_result

    # ==========================================================================
    # Stage Implementations
    # ==========================================================================

    def _stage_initialization(self) -> Dict:
        """Initialize workspace and validate configuration."""
        outputs = {}
        warnings = []

        # Validate VOC name
        voc_name = self.config.voc_name.lower().strip()
        if not voc_name:
            raise ValueError("VOC name cannot be empty")

        # Check for VOC in database
        try:
            from gecko_web.chemdata import get_compound
            compound = get_compound(voc_name)
            if compound:
                outputs['compound_info'] = {
                    'name': compound.name,
                    'smiles': compound.smiles,
                    'gecko_formula': compound.gecko_formula,
                    'molecular_weight': compound.molecular_weight,
                    'category': compound.category
                }
                self._log(f"Found compound in database: {compound.name}")
            else:
                warnings.append(f"VOC '{voc_name}' not found in compound database")
        except ImportError:
            warnings.append("Compound database not available")

        # Write configuration
        config_file = self.output.reports_dir / "workflow_config.json"
        config_dict = asdict(self.config)
        # Convert Path objects to strings for JSON serialization
        for key, value in config_dict.items():
            if isinstance(value, Path):
                config_dict[key] = str(value)
        with open(config_file, 'w') as f:
            json.dump(config_dict, f, indent=2, default=str)

        outputs['config_file'] = str(config_file)

        return {
            'message': 'Initialization complete',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_mechanism_generation(self) -> Dict:
        """Run GECKO-A mechanism generator."""
        outputs = {}
        warnings = []

        # Get VOC formula
        voc_formula = self.config.voc_formula
        if not voc_formula:
            try:
                from gecko_web.chemdata import get_gecko_formula
                voc_formula = get_gecko_formula(self.config.voc_name)
            except (ImportError, Exception):
                pass

        if not voc_formula:
            # Fall back to VOC mapping in main
            try:
                from gecko_web.main import get_gecko_input
                voc_formula = get_gecko_input(self.config.voc_name)
            except (ImportError, Exception):
                voc_formula = self.config.voc_name

        self._log(f"Using GECKO formula: {voc_formula}")

        # Check environment - support multiple possible locations
        gecko_root = None
        possible_paths = [
            Path(os.getenv("GECKO_SOURCE_DIR", "")),
            Path("/app/gecko_source"),
            Path.home() / "gecko_source",
            Path("/usr/local/gecko"),
        ]

        for path in possible_paths:
            if path.exists() and (path / "RUN").exists():
                gecko_root = path
                break

        if gecko_root is None:
            # Check if we have a library cache that can be used
            try:
                from gecko_web.main import DATA_DIR
                library_path = DATA_DIR / "library" / self.config.voc_name.lower().strip()
                if library_path.exists():
                    self._log(f"Using cached mechanism from {library_path}")
                    for item in library_path.iterdir():
                        if item.is_file():
                            shutil.copy2(item, self.output.mechanism_dir)
                    outputs['from_cache'] = True
                    outputs['mechanism_dir'] = str(self.output.mechanism_dir)

                    # Count species from cache
                    dict_file = self.output.mechanism_dir / "dictionary.out"
                    if dict_file.exists():
                        with open(dict_file) as f:
                            species_count = sum(1 for line in f if line.strip()
                                               and not line.startswith('!')
                                               and 'number of record' not in line.lower())
                        outputs['species_count'] = species_count

                    return {
                        'message': f'Loaded cached mechanism with {outputs.get("species_count", "?")} species',
                        'outputs': outputs,
                        'warnings': warnings
                    }
            except ImportError:
                pass

            raise EnvironmentError(
                f"GECKO-A not found. Searched paths: {[str(p) for p in possible_paths]}. "
                "Set GECKO_SOURCE_DIR environment variable or ensure GECKO-A is installed."
            )

        # Import generator functions from main
        try:
            from gecko_web.main import (
                run_gecko_generator,
                workspace_manager,
                job_manager,
                ScenarioParams
            )

            # Use the workflow_id as the job_id for tracking
            # This avoids creating a separate temp job that can cause confusion
            job_id = self.workflow_id

            # Create scenario params
            scenario = ScenarioParams(
                vapor_pressure_threshold=self.config.vapor_pressure_threshold,
                max_generations=self.config.max_generations,
                temperature_k=self.config.temperature_k
            )

            # Run generator with workflow_id for logging
            self._log(f"Running GECKO-A generator...")
            success = run_gecko_generator(
                job_id,
                self.config.voc_name,
                self.output.mechanism_dir,
                scenario
            )

            if success is None or success is False:
                raise RuntimeError("GECKO-A generator returned failure or no result")

            # Count output files
            mech_files = list(self.output.mechanism_dir.glob("*"))
            outputs['files_generated'] = len(mech_files)
            outputs['mechanism_dir'] = str(self.output.mechanism_dir)

            # Check for key files
            dict_file = self.output.mechanism_dir / "dictionary.out"
            if dict_file.exists():
                # Count species
                with open(dict_file) as f:
                    species_count = sum(1 for line in f if line.strip()
                                       and not line.startswith('!')
                                       and 'number of record' not in line.lower())
                outputs['species_count'] = species_count
                self._log(f"Generated mechanism with {species_count} species")
            else:
                warnings.append("dictionary.out not found in mechanism output")

            # Copy mechanism files to exports
            for mec_file in self.output.mechanism_dir.glob("*.mec"):
                shutil.copy2(mec_file, self.output.exports_dir)

        except ImportError as e:
            self._log(f"Import error: {e}", level="error")
            warnings.append(f"Could not import GECKO-A functions: {e}")
            raise RuntimeError(f"GECKO-A generator not available: {e}")

        return {
            'message': f'Mechanism generated with {outputs.get("species_count", "?")} species',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_mass_balance(self) -> Dict:
        """Verify mass balance of generated mechanism."""
        outputs = {}
        warnings = []

        try:
            from gecko_web.mass_balance import (
                verify_mass_balance,
                generate_mass_balance_json,
                save_mass_balance_report
            )

            # Run mass balance verification
            report = verify_mass_balance(
                str(self.output.mechanism_dir),
                self.config.voc_name,
                verbose=self.config.verbose
            )

            # Save reports
            json_path = self.output.analysis_dir / "mass_balance.json"
            txt_path = self.output.analysis_dir / "mass_balance.txt"

            save_mass_balance_report(report, str(json_path), format='json')
            save_mass_balance_report(report, str(txt_path), format='txt')

            outputs['mass_balance_json'] = str(json_path)
            outputs['mass_balance_txt'] = str(txt_path)
            outputs['is_valid'] = report.is_valid
            outputs['balance_percentage'] = (
                report.balanced_reactions / report.total_reactions * 100
                if report.total_reactions > 0 else 0
            )
            outputs['critical_errors'] = report.critical_errors

            if report.critical_errors > 0:
                warnings.append(
                    f"Mass balance check found {report.critical_errors} critical errors"
                )

            if not report.is_valid:
                warnings.extend(report.recommendations)

            self._log(f"Mass balance: {outputs['balance_percentage']:.1f}% balanced")

        except ImportError:
            warnings.append("Mass balance module not available")
            outputs['is_valid'] = None

        return {
            'message': f"Mass balance check: {outputs.get('balance_percentage', 'N/A'):.1f}% balanced",
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_box_model(self) -> Dict:
        """Run BOXMODEL4GECKO simulation."""
        outputs = {}
        warnings = []

        boxmodel_root = Path(os.getenv("BOXMODEL_SOURCE_DIR", "/app/boxmodel_source"))
        if not boxmodel_root.exists():
            warnings.append(f"BOXMODEL not found at {boxmodel_root}")
            return {
                'message': 'Box model skipped (not available)',
                'outputs': outputs,
                'warnings': warnings
            }

        try:
            from gecko_web.main import (
                run_box_model_with_mechanism_dir,
                job_manager,
                ScenarioParams
            )

            # Create scenario
            scenario = ScenarioParams(
                temperature_k=self.config.temperature_k,
                rh_percent=self.config.rh_percent,
                latitude_degrees=self.config.latitude_degrees,
                initial_o3_ppb=self.config.initial_o3_ppb,
                initial_nox_ppb=self.config.initial_nox_ppb,
                seed_aerosol_ug_m3=self.config.seed_aerosol_ug_m3,
                dilution_rate_s1=self.config.dilution_rate_s1
            )

            # Create temp job for logging
            temp_job = job_manager.create_job(
                "boxmodel",
                self.config.voc_name,
                {'workflow_id': self.workflow_id}
            )

            # Run box model with explicit mechanism directory
            # The mechanism files are in mechanism_dir, results go to simulation_dir
            success = run_box_model_with_mechanism_dir(
                temp_job.id,
                self.config.voc_name,
                mechanism_dir=self.output.mechanism_dir,
                output_dir=self.output.simulation_dir,
                scenario=scenario
            )

            if not success:
                raise RuntimeError("Box model returned failure")

            # Count output files
            sim_files = list(self.output.simulation_dir.glob("*"))
            outputs['files_generated'] = len(sim_files)
            outputs['simulation_dir'] = str(self.output.simulation_dir)

            # Check for netCDF output
            nc_files = list(self.output.simulation_dir.glob("*.nc"))
            if nc_files:
                outputs['netcdf_file'] = str(nc_files[0])

        except ImportError as e:
            warnings.append(f"Could not import box model functions: {e}")
            raise RuntimeError("Box model not available")

        return {
            'message': f'Box model simulation complete',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_post_processing(self) -> Dict:
        """Run post-processing analysis."""
        outputs = {}
        warnings = []

        try:
            from gecko_web import postprocessing

            # Determine which directory has the data
            if list(self.output.simulation_dir.glob("*.nc")):
                data_dir = self.output.simulation_dir
            else:
                data_dir = self.output.mechanism_dir

            # Run post-processing
            result = postprocessing.run_postprocessing(
                str(data_dir),
                self.config.voc_name,
                temperature_k=self.config.temperature_k,
                seed_mass=self.config.seed_aerosol_ug_m3
            )

            if result['status'] == 'error':
                warnings.append(f"Post-processing error: {result.get('error')}")
            else:
                # Move generated files to appropriate directories
                for key, path in result.get('outputs', {}).items():
                    if path and os.path.exists(path):
                        src = Path(path)
                        if src.suffix == '.csv':
                            dst = self.output.analysis_dir / src.name
                        elif src.suffix == '.json':
                            dst = self.output.analysis_dir / src.name
                        elif src.suffix == '.png':
                            dst = self.output.plots_dir / src.name
                        else:
                            continue

                        shutil.copy2(src, dst)
                        outputs[key] = str(dst)

                # Add summary data
                if result.get('summary'):
                    outputs['summary'] = result['summary']

            warnings.extend(result.get('warnings', []))

        except ImportError:
            warnings.append("Post-processing module not available")
        except Exception as e:
            warnings.append(f"Post-processing failed: {e}")

        return {
            'message': 'Post-processing complete',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_diagram_generation(self) -> Dict:
        """Generate mechanism diagrams."""
        outputs = {}
        warnings = []

        try:
            from gecko_web import reaction_tree, mechanism_diagram, pathway_visualizer, enhanced_visualizer

            # Parse reaction tree
            tree_data = reaction_tree.parse_reaction_tree(
                str(self.output.mechanism_dir),
                self.config.voc_name,
                max_depth=self.config.diagram_max_depth,
                max_nodes=self.config.max_diagram_nodes
            )

            # Save reaction tree JSON
            tree_file = self.output.diagrams_dir / "reaction_tree.json"
            with open(tree_file, 'w') as f:
                json.dump(tree_data, f, indent=2)
            outputs['reaction_tree'] = str(tree_file)

            # Generate pathway diagram using Enhanced Visualizer (Color-Coded)
            if tree_data.get('nodes'):
                try:
                    # Use enhanced visualizer for publication-quality colorful diagrams
                    visualizer = enhanced_visualizer.EnhancedPathwayVisualizer()
                    visualizer.from_tree_data(tree_data)
                    
                    # Generate PNG
                    visualizer.generate_diagram(
                        output_path=str(self.output.diagrams_dir / "pathway_diagram.png"),
                        title=f"Oxidation Mechanism: {self.config.voc_name.upper()}",
                        layout='hierarchical',
                        show_branching_ratios=True,
                        color_by_reaction_type=True,
                        format='png'
                    )
                    outputs['pathway_png'] = str(self.output.diagrams_dir / "pathway_diagram.png")
                    
                    # Generate SVG
                    visualizer.generate_diagram(
                         output_path=str(self.output.diagrams_dir / "pathway_diagram.svg"),
                         title=f"Oxidation Mechanism: {self.config.voc_name.upper()}",
                         layout='hierarchical',
                         show_branching_ratios=True,
                         color_by_reaction_type=True,
                         format='svg'
                    )
                    outputs['pathway_svg'] = str(self.output.diagrams_dir / "pathway_diagram.svg")

                except Exception as e:
                    logger.warning(f"Enhanced pathway diagram generation failed ({e}), falling back to legacy visualizer...")
                    try:
                        pathway_visualizer.generate_pathway_diagram(
                            tree_data=tree_data,
                            output_path=str(self.output.diagrams_dir / "pathway_diagram"),
                            title=f"Oxidation Mechanism: {self.config.voc_name.upper()}",
                            format='png',
                            mol_size=(180, 140)
                        )
                        outputs['pathway_png'] = str(self.output.diagrams_dir / "pathway_diagram.png")
                    except Exception as e2:
                        warnings.append(f"Pathway diagram generation completely failed: {e2}")

            # Generate mechanism exports
            if self.config.generate_all_formats:
                try:
                    mechanism_diagram.generate_all_outputs(
                        output_dir=str(self.output.exports_dir),
                        voc_name=self.config.voc_name,
                        generate_kpp=True,
                        generate_mcm=True,
                        generate_facsimile=True
                    )
                    outputs['exports_generated'] = True
                except Exception as e:
                    warnings.append(f"Export generation failed: {e}")

        except ImportError as e:
            warnings.append(f"Diagram modules not available: {e}")
        except Exception as e:
            warnings.append(f"Diagram generation failed: {e}")

        return {
            'message': 'Diagrams generated',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_report_generation(self) -> Dict:
        """Generate summary reports."""
        outputs = {}
        warnings = []

        # Generate text summary
        summary_lines = [
            f"GECKO-A Workflow Report",
            f"=" * 60,
            f"VOC: {self.config.voc_name}",
            f"Workflow ID: {self.workflow_id}",
            f"Generated: {datetime.now().isoformat()}",
            "",
            "Configuration:",
            f"  Max Generations: {self.config.max_generations}",
            f"  VP Threshold: {self.config.vapor_pressure_threshold}",
            f"  Temperature: {self.config.temperature_k} K",
            "",
            "Stage Results:",
        ]

        for stage in self.stages:
            status_symbol = "✓" if stage.status == WorkflowStatus.COMPLETED else "✗"
            summary_lines.append(
                f"  {status_symbol} {stage.stage.value}: {stage.status.value} "
                f"({stage.duration_seconds:.2f}s)"
            )

        # Add analysis summary if available
        analysis_summary = self.output.analysis_dir / "postprocessing_summary.json"
        if analysis_summary.exists():
            try:
                with open(analysis_summary) as f:
                    data = json.load(f)
                summary_lines.extend([
                    "",
                    "Analysis Results:",
                    f"  Total SOA Mass: {data.get('total_soa_mass_ug_m3', 'N/A'):.3f} µg/m³",
                    f"  SOA Yield: {data.get('soa_yield', 'N/A'):.4f}",
                    f"  C* (median): {data.get('median_log_c_star', 'N/A'):.2f}"
                ])
            except Exception:
                pass

        summary_text = "\n".join(summary_lines)

        # Save text summary
        summary_file = self.output.reports_dir / "summary.txt"
        with open(summary_file, 'w') as f:
            f.write(summary_text)
        outputs['summary_txt'] = str(summary_file)

        # Generate PDF report if enabled
        if self.config.generate_pdf_report:
            try:
                pdf_path = self._generate_pdf_report()
                if pdf_path:
                    outputs['pdf_report'] = str(pdf_path)
            except Exception as e:
                warnings.append(f"PDF generation failed: {e}")

        return {
            'message': 'Reports generated',
            'outputs': outputs,
            'warnings': warnings
        }

    def _stage_finalization(self) -> Dict:
        """Finalize workflow and cleanup."""
        outputs = {}
        warnings = []

        # List all outputs
        all_outputs = self.output.list_outputs()
        outputs['output_summary'] = all_outputs

        total_files = sum(len(files) for files in all_outputs.values())
        outputs['total_files'] = total_files

        # Calculate total size
        total_size = 0
        for subdir in self.output.SUBDIRS:
            subdir_path = self.output.root / subdir
            if subdir_path.exists():
                for f in subdir_path.rglob("*"):
                    if f.is_file():
                        total_size += f.stat().st_size

        outputs['total_size_mb'] = total_size / (1024 * 1024)

        # Cleanup intermediate files if not preserving
        if not self.config.preserve_intermediate:
            # Could remove temp files here
            pass

        return {
            'message': f'Workflow finalized: {total_files} files, {outputs["total_size_mb"]:.2f} MB',
            'outputs': outputs,
            'warnings': warnings
        }

    def _build_summary(self) -> Dict:
        """Build workflow summary."""
        summary = {
            'voc_name': self.config.voc_name,
            'total_stages': len(self.stages),
            'completed_stages': sum(
                1 for s in self.stages
                if s.status in (WorkflowStatus.COMPLETED, WorkflowStatus.WARNING)
            ),
            'failed_stages': sum(
                1 for s in self.stages if s.status == WorkflowStatus.FAILED
            ),
            'total_warnings': sum(len(s.warnings) for s in self.stages),
            'outputs': self.output.list_outputs()
        }

        # Add key metrics if available
        analysis_summary = self.output.analysis_dir / "postprocessing_summary.json"
        if analysis_summary.exists():
            try:
                with open(analysis_summary) as f:
                    data = json.load(f)
                summary['analysis'] = data
            except Exception:
                pass

        return summary

    def _generate_pdf_report(self) -> Optional[Path]:
        """Generate PDF report with embedded plots."""
        try:
            from reportlab.lib.pagesizes import letter
            from reportlab.platypus import (
                SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
            )
            from reportlab.lib.styles import getSampleStyleSheet
            from reportlab.lib import colors
        except ImportError:
            self._log("ReportLab not available - skipping PDF generation")
            return None

        pdf_path = self.output.reports_dir / "workflow_report.pdf"

        doc = SimpleDocTemplate(str(pdf_path), pagesize=letter)
        styles = getSampleStyleSheet()
        elements = []

        # Title
        elements.append(Paragraph(
            f"GECKO-A Workflow Report: {self.config.voc_name}",
            styles['Title']
        ))
        elements.append(Spacer(1, 20))

        # Metadata
        elements.append(Paragraph("Workflow Information", styles['Heading2']))
        metadata = [
            ["Workflow ID", self.workflow_id],
            ["VOC Name", self.config.voc_name],
            ["Date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
            ["Max Generations", str(self.config.max_generations)],
            ["Temperature", f"{self.config.temperature_k} K"],
        ]
        t = Table(metadata, colWidths=[150, 300])
        t.setStyle(TableStyle([
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
        ]))
        elements.append(t)
        elements.append(Spacer(1, 20))

        # Stage results
        elements.append(Paragraph("Stage Results", styles['Heading2']))
        stage_data = [["Stage", "Status", "Duration (s)"]]
        for stage in self.stages:
            stage_data.append([
                stage.stage.value,
                stage.status.value,
                f"{stage.duration_seconds:.2f}"
            ])
        t = Table(stage_data, colWidths=[200, 100, 100])
        t.setStyle(TableStyle([
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
        ]))
        elements.append(t)
        elements.append(Spacer(1, 20))

        # Add plots
        elements.append(Paragraph("Analysis Plots", styles['Heading2']))

        plot_files = list(self.output.plots_dir.glob("*.png"))
        for plot_file in plot_files[:6]:  # Limit to 6 plots
            try:
                img = Image(str(plot_file), width=400, height=300)
                elements.append(img)
                elements.append(Paragraph(plot_file.stem.replace("_", " ").title(), styles['Caption']))
                elements.append(Spacer(1, 10))
            except Exception:
                pass

        doc.build(elements)
        return pdf_path


# ==============================================================================
# Convenience Function
# ==============================================================================

def run_combined_workflow(
    voc_name: str,
    output_dir: str,
    run_box_model: bool = True,
    max_generations: int = 2,
    temperature_k: float = 298.0,
    verbose: bool = False,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> WorkflowResult:
    """
    Run the combined GECKO-A workflow with default settings.

    Args:
        voc_name: Name of the VOC (e.g., "alpha-pinene")
        output_dir: Base directory for outputs
        run_box_model: Whether to run box model simulation
        max_generations: Maximum reaction generations
        temperature_k: Simulation temperature
        verbose: Print progress to console
        progress_callback: Optional progress callback

    Returns:
        WorkflowResult with execution details
    """
    config = WorkflowConfig(
        voc_name=voc_name,
        run_box_model=run_box_model,
        max_generations=max_generations,
        temperature_k=temperature_k,
        verbose=verbose
    )

    workflow = CombinedWorkflow(
        config=config,
        output_dir=Path(output_dir),
        progress_callback=progress_callback
    )

    return workflow.run()


# ==============================================================================
# CLI Entry Point
# ==============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Run GECKO-A combined workflow')
    parser.add_argument('voc_name', help='VOC name (e.g., alpha-pinene)')
    parser.add_argument('--output', '-o', default='./output', help='Output directory')
    parser.add_argument('--no-boxmodel', action='store_true', help='Skip box model')
    parser.add_argument('--generations', '-g', type=int, default=2, help='Max generations')
    parser.add_argument('--temperature', '-T', type=float, default=298.0, help='Temperature (K)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    def progress_cb(msg, frac):
        print(f"[{frac*100:.0f}%] {msg}")

    result = run_combined_workflow(
        voc_name=args.voc_name,
        output_dir=args.output,
        run_box_model=not args.no_boxmodel,
        max_generations=args.generations,
        temperature_k=args.temperature,
        verbose=args.verbose,
        progress_callback=progress_cb if args.verbose else None
    )

    print(f"\nWorkflow completed: {result.status.value}")
    print(f"Output directory: {result.output_directory}")
    print(f"Total duration: {result.total_duration_seconds:.2f}s")
