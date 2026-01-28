"""
GECKO-A Post-Processing Module - Scientifically Rigorous Implementation

This module performs aerosol partitioning calculations and generates analysis plots.

KEY SCIENTIFIC IMPROVEMENTS:
1. Dynamic C_OA calculation: Iteratively solves the partitioning equations
   at each timestep instead of using a static value
2. Smart Fortran Splitter: Correctly parses merged Fortran output columns
3. Real Vapor Pressure: Reads pvap.nan.dat for actual calculated vapor pressures
   instead of crude MW-based correlations

References:
- Pankow (1994): Absorptive partitioning theory
- Odum et al. (1996): SOA yield parameterization
- Donahue et al. (2006): Volatility Basis Set

Author: Deeksha Sharma
"""

import os
import re
import glob
import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Physical constants
R_GAS = 8.314  # J/(mol·K)
R_ATM = 8.206e-5  # m³·atm/(mol·K)
MW_AIR = 29.0  # g/mol


# ==============================================================================
# Smart Fortran Output Parser
# ==============================================================================

class FortranOutputParser:
    """
    Handles parsing of Fortran fixed-width output files.

    Fortran's fixed-width formatting causes numbers to merge when they overflow
    their allocated column width. For example:
        Expected: "  120.45   3.2E-5"
        Actual:   "  120.453.2E-5" (columns merged)

    This parser intelligently splits merged numbers by detecting the pattern
    of scientific notation boundaries.
    """

    # Pattern to match scientific notation numbers with LIMITED exponent digits (1-2)
    # This prevents greedy matching of merged numbers
    SCIENTIFIC_PATTERN = re.compile(
        r'[+-]?\d+\.?\d*[EeDd][+-]?\d{1,2}(?!\d)'
    )

    # Pattern for standard floats
    FLOAT_PATTERN = re.compile(
        r'[+-]?\d+\.\d+'
    )

    @classmethod
    def smart_split(cls, text: str) -> List[str]:
        """
        Intelligently split a string that may contain merged Fortran numbers.

        Args:
            text: Input string potentially with merged numbers

        Returns:
            List of individual number strings
        """
        text = text.strip()
        if not text:
            return []

        # First try: standard whitespace split
        parts = text.split()

        # Check if any part needs further splitting
        result = []
        for part in parts:
            split_numbers = cls._split_merged_numbers(part)
            result.extend(split_numbers)

        return result

    @classmethod
    def _split_merged_numbers(cls, text: str) -> List[str]:
        """
        Split a string that may contain merged numbers.

        Handles cases like:
        - "1.00E+012.00E+01" -> ["1.00E+01", "2.00E+01"]
        - "120.453.2E-5" -> ["120.45", "3.2E-5"]
        - "1.5E-102.3E-10" -> ["1.5E-10", "2.3E-10"]
        """
        if not text:
            return []

        # Special handling for merged scientific notation
        # Look for patterns like number.numberE or E+digit.digit
        numbers = cls._split_scientific_notation(text)
        if numbers and len(numbers) > 1:
            return numbers

        # Try to extract all scientific notation numbers
        sci_matches = list(cls.SCIENTIFIC_PATTERN.finditer(text))

        if sci_matches:
            numbers = []
            last_end = 0

            for match in sci_matches:
                # Check if there's a number before this scientific number
                prefix = text[last_end:match.start()]
                if prefix:
                    # There might be a regular float before the sci notation
                    prefix_numbers = cls._extract_regular_floats(prefix)
                    numbers.extend(prefix_numbers)

                numbers.append(match.group())
                last_end = match.end()

            # Check for trailing numbers
            suffix = text[last_end:]
            if suffix:
                suffix_numbers = cls._extract_regular_floats(suffix)
                numbers.extend(suffix_numbers)

            if numbers:
                return numbers

        # No scientific notation, try regular floats
        float_numbers = cls._extract_regular_floats(text)
        if float_numbers:
            return float_numbers

        # Last resort: return original if it looks like a single number
        try:
            float(text)
            return [text]
        except ValueError:
            return [text]  # Return as-is for non-numeric tokens

    @classmethod
    def _split_scientific_notation(cls, text: str) -> List[str]:
        """
        Split merged scientific notation numbers by looking for E/D patterns.

        Fortran scientific notation typically has 1-2 digit exponents.
        When numbers merge, we see patterns like:
        - 1.00E+012.00E+01 (should be 1.00E+01, 2.00E+01)
        - 1.5E-102.3E-10 (should be 1.5E-10, 2.3E-10)
        """
        if not text:
            return []

        # Look for E/D with sign and 2+ digits (indicates merged numbers)
        # Pattern: digit after E+/-NN where there are more than 2 digits
        merge_pattern = re.compile(r'([EeDd][+-]?)(\d{3,})')

        match = merge_pattern.search(text)
        if match:
            exp_sign = match.group(1)
            digits = match.group(2)

            # Fortran exponents are typically 1-2 digits
            # Split at position 2 (e.g., "012" -> "01" and "2...")
            if len(digits) >= 3:
                # The first 2 digits belong to first number's exponent
                exp_digits = digits[:2]
                remainder = digits[2:]

                # Reconstruct first number
                first_num = text[:match.start()] + exp_sign + exp_digits

                # Parse remainder recursively
                remaining_text = remainder + text[match.end():]
                other_nums = cls._split_merged_numbers(remaining_text)

                return [first_num] + other_nums

        return []

    @classmethod
    def _extract_regular_floats(cls, text: str) -> List[str]:
        """Extract regular floating point numbers from text."""
        matches = cls.FLOAT_PATTERN.findall(text)
        return matches if matches else []

    @classmethod
    def parse_value(cls, text: str, default: float = 0.0) -> float:
        """
        Parse a single value from potentially merged text.

        Args:
            text: Input text
            default: Default value if parsing fails

        Returns:
            Parsed float value
        """
        numbers = cls.smart_split(text)
        if numbers:
            try:
                return float(numbers[0])
            except ValueError:
                pass
        return default


# ==============================================================================
# Vapor Pressure Reader
# ==============================================================================

class VaporPressureReader:
    """
    Reads vapor pressure data from GECKO-A output files.

    GECKO-A calculates vapor pressures using sophisticated SARs (Structure-Activity
    Relationships) like Nannoolal or SIMPOL-1. These are stored in pvap.nan.dat.

    Using these real values instead of crude MW correlations improves accuracy
    by orders of magnitude.
    """

    @staticmethod
    def read_pvap_file(output_dir: str) -> Dict[str, float]:
        """
        Read vapor pressure data from pvap.nan.dat.

        File format (space-separated):
        CODE  Pvap_atm  Hvap_kJ_mol  [other fields...]

        The file contains actual Pvap values in atm (e.g., 2.339E-06), not log10 values.
        This function converts to log10(Pvap) for internal use.

        Returns:
            Dict mapping species code to log10(Pvap) in atm
        """
        pvap_data = {}

        pvap_files = [
            os.path.join(output_dir, "pvap.nan.dat"),
            os.path.join(output_dir, "pvap.sim.dat"),
            os.path.join(output_dir, "pvap.jrm.dat"),
        ]

        for pvap_file in pvap_files:
            if os.path.exists(pvap_file):
                logger.info(f"Reading vapor pressures from {pvap_file}")
                try:
                    with open(pvap_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            if not line or line.startswith('!') or line.startswith('#'):
                                continue
                            # Skip TREF line
                            if line.startswith('TREF'):
                                continue

                            parts = FortranOutputParser.smart_split(line)
                            if len(parts) >= 2:
                                code = parts[0]
                                try:
                                    pvap_atm = float(parts[1])
                                    # Convert to log10(Pvap) - file contains actual Pvap, not log
                                    if pvap_atm > 0:
                                        log_pvap = np.log10(pvap_atm)
                                    else:
                                        log_pvap = -30  # Essentially zero vapor pressure
                                    pvap_data[code] = log_pvap
                                    # Also store without 'G' prefix
                                    if code.startswith('G'):
                                        pvap_data[code[1:]] = log_pvap
                                except ValueError:
                                    continue
                    break  # Use first available file
                except Exception as e:
                    logger.warning(f"Error reading {pvap_file}: {e}")

        logger.info(f"Loaded vapor pressures for {len(pvap_data)} species")
        return pvap_data

    @staticmethod
    def estimate_from_mw(mw: float, n_oxygen: int = 0, n_nitrogen: int = 0) -> float:
        """
        Fallback: Estimate log10(Pvap) from molecular weight and composition.

        This uses a simple correlation based on functional group contributions.
        Much less accurate than SAR methods but serves as a fallback.

        Based on: Modified Clausius-Clapeyron with functional group corrections.
        """
        if mw <= 0:
            return 0.0

        # Base correlation (hydrocarbons)
        # log10(P_L^0) ~ A - B * MW
        # Fitted to alkanes: ~6.5 - 0.04 * MW

        log_pvap = 6.5 - 0.04 * mw

        # Oxygen correction: each O reduces volatility ~1 order of magnitude
        # This accounts for H-bonding and polarity
        log_pvap -= 0.8 * n_oxygen

        # Nitrogen correction (smaller effect)
        log_pvap -= 0.3 * n_nitrogen

        return log_pvap


# ==============================================================================
# Dictionary Parser with Smart Splitting
# ==============================================================================

@dataclass
class SpeciesData:
    """Data structure for a chemical species."""
    code: str
    formula: str
    mw: float
    n_carbon: int
    n_hydrogen: int
    n_nitrogen: int
    n_oxygen: int
    n_sulfur: int = 0
    log_pvap: Optional[float] = None
    c_star: Optional[float] = None


def parse_dictionary_robust(dictionary_path: str) -> Dict[str, SpeciesData]:
    """
    Parse dictionary.out with robust handling of Fortran formatting issues.

    The dictionary file has variable formatting that can cause column merging.
    This parser handles those cases gracefully.
    """
    species_dict = {}

    if not os.path.exists(dictionary_path):
        logger.warning(f"Dictionary not found: {dictionary_path}")
        return species_dict

    with open(dictionary_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('!') or 'number of record' in line.lower():
                continue

            try:
                # Use smart split for the entire line
                parts = FortranOutputParser.smart_split(line)
                if len(parts) < 4:
                    continue

                code = parts[0]
                formula = parts[1]

                # Find the molecular weight (first float with decimal that looks like MW)
                # MW values are typically 10-1000, class codes like "1." or "2." are too small
                mw = 0.0
                mw_idx = -1
                for idx, p in enumerate(parts[2:], start=2):
                    try:
                        val = float(p)
                        # Must have decimal AND be a reasonable MW (>10 Da)
                        if ('.' in p or 'E' in p.upper()) and val >= 10.0:
                            mw = val
                            mw_idx = idx
                            break
                    except ValueError:
                        continue

                # Element counts follow the MW
                n_carbon, n_hydrogen, n_nitrogen, n_oxygen, n_sulfur = 0, 0, 0, 0, 0

                if mw_idx >= 0:
                    # After MW, there's typically a radical/generation flag (0 or 1), then element counts
                    # Dictionary format: code, formula, classCode, MW, radicalFlag, nC, nH, nN, nO, nS, ...
                    elem_start = mw_idx + 1

                    # Skip phase flag if present (G, P, A)
                    if elem_start < len(parts) and parts[elem_start] in ('G', 'P', 'A'):
                        elem_start += 1

                    # Skip radical/generation flag (0 or 1) - this is NOT the carbon count!
                    # The flag indicates if the species is a radical (1) or closed-shell (0)
                    if elem_start < len(parts):
                        try:
                            flag = int(parts[elem_start])
                            if flag in (0, 1):
                                elem_start += 1
                        except ValueError:
                            pass

                    # Extract element counts (order: nC, nH, nN, nO, nS typically)
                    elem_values = []
                    for p in parts[elem_start:elem_start + 6]:
                        try:
                            elem_values.append(int(float(p)))
                        except ValueError:
                            break

                    if len(elem_values) >= 1:
                        n_carbon = elem_values[0]
                    if len(elem_values) >= 2:
                        n_hydrogen = elem_values[1]
                    if len(elem_values) >= 3:
                        n_nitrogen = elem_values[2]
                    if len(elem_values) >= 4:
                        n_oxygen = elem_values[3]
                    if len(elem_values) >= 5:
                        n_sulfur = elem_values[4]

                # Validate MW
                if mw < 1.0 or mw > 10000.0:
                    # Try to estimate from formula
                    mw = _estimate_mw_from_composition(n_carbon, n_hydrogen, n_nitrogen, n_oxygen, n_sulfur)

                species_dict[code] = SpeciesData(
                    code=code,
                    formula=formula,
                    mw=mw,
                    n_carbon=n_carbon,
                    n_hydrogen=n_hydrogen,
                    n_nitrogen=n_nitrogen,
                    n_oxygen=n_oxygen,
                    n_sulfur=n_sulfur
                )

                # Also store with 'G' prefix
                species_dict[f'G{code}'] = species_dict[code]

            except Exception as e:
                logger.debug(f"Failed to parse line {line_num}: {e}")
                continue

    logger.info(f"Parsed {len(species_dict)} species from dictionary")
    return species_dict


def _estimate_mw_from_composition(n_c: int, n_h: int, n_n: int, n_o: int, n_s: int = 0) -> float:
    """Estimate molecular weight from elemental composition."""
    return n_c * 12.011 + n_h * 1.008 + n_n * 14.007 + n_o * 15.999 + n_s * 32.065


# ==============================================================================
# Dynamic Gas-Particle Partitioning Calculator
# ==============================================================================

class PartitioningCalculator:
    """
    Implements dynamic gas-particle partitioning calculations.

    The key equation (Pankow, 1994):
        F_p,i = C_OA / (C_OA + C*_i)

    where C_OA is the TOTAL organic aerosol mass concentration, which itself
    depends on the partitioning:
        C_OA = sum_i(C_i,total * F_p,i)

    This creates a coupled system that requires iterative solution.
    """

    GAS_CONSTANT = 8.314  # J/(mol·K)
    ACTIVITY_COEFFICIENT = 1.0  # Assume ideal mixing (Raoult's law)

    def __init__(self, temperature_k: float = 298.0, seed_mass: float = 0.0):
        """
        Initialize calculator.

        Args:
            temperature_k: Temperature in Kelvin
            seed_mass: Initial seed aerosol mass (ug/m3) - inert, adds to C_OA
        """
        self.temperature_k = temperature_k
        self.seed_mass = seed_mass

    def calculate_c_star(self, mw: float, log_pvap_atm: float) -> float:
        """
        Calculate saturation concentration C* (ug/m3).

        C* = (MW * Pvap_Pa) / (R * T) * 10^6

        where:
        - MW is in g/mol
        - Pvap_Pa is vapor pressure in Pa (converted from atm: 1 atm = 101325 Pa)
        - R = 8.314 J/(mol·K) = 8.314 Pa·m³/(mol·K)
        - T is in K
        - Result is in g/m³, multiplied by 10^6 to get µg/m³

        References:
        - Pankow (1994): Atmos. Environ. 28, 185-188
        - Donahue et al. (2006): Environ. Sci. Technol. 40, 2635-2643
        """
        if mw <= 0 or log_pvap_atm > 10:  # Sanity check
            return 1e10  # Very volatile

        pvap_atm = 10 ** log_pvap_atm
        pvap_pa = pvap_atm * 101325  # Convert atm to Pa

        # C* in ug/m3
        # C* = MW * Pvap / (R * T) gives g/m³, multiply by 10^6 for µg/m³
        c_star = (mw * pvap_pa) / (self.GAS_CONSTANT * self.temperature_k) * 1e6

        return max(c_star, 1e-10)  # Prevent division by zero

    def calculate_partitioning_static(self, species_data: List[Dict],
                                      c_oa_assumed: float = 10.0) -> pd.DataFrame:
        """
        Calculate partitioning with static C_OA (for comparison/legacy).

        This is the WRONG approach but included for comparison.
        """
        results = []
        for sp in species_data:
            if sp['c_star'] > 0:
                f_p = c_oa_assumed / (c_oa_assumed + sp['c_star'])
            else:
                f_p = 1.0

            results.append({
                'code': sp['code'],
                'mw': sp['mw'],
                'c_star': sp['c_star'],
                'c_total': sp['concentration'],
                'f_p': f_p,
                'c_particle': sp['concentration'] * f_p,
                'c_gas': sp['concentration'] * (1 - f_p),
                'c_oa_used': c_oa_assumed
            })

        return pd.DataFrame(results)

    def calculate_partitioning_dynamic(self, species_data: List[Dict],
                                       max_iterations: int = 100,
                                       tolerance: float = 1e-6) -> pd.DataFrame:
        """
        Calculate partitioning with dynamic C_OA (CORRECT approach).

        Uses iterative solution:
        1. Start with initial C_OA estimate (seed mass + estimate from volatilities)
        2. Calculate F_p for all species using current C_OA
        3. Calculate new C_OA = seed + sum(C_i * F_p,i)
        4. Repeat until convergence

        Args:
            species_data: List of dicts with 'code', 'mw', 'c_star', 'concentration'
            max_iterations: Maximum iteration count
            tolerance: Convergence tolerance (relative change in C_OA)

        Returns:
            DataFrame with partitioning results
        """
        n_species = len(species_data)
        if n_species == 0:
            return pd.DataFrame()

        # Extract arrays for vectorized operations
        c_stars = np.array([sp['c_star'] for sp in species_data])
        c_totals = np.array([sp['concentration'] for sp in species_data])

        # Handle zero/invalid concentrations
        c_totals = np.where(c_totals > 0, c_totals, 0.0)
        c_stars = np.where(c_stars > 0, c_stars, 1e-10)

        # Initial C_OA estimate
        # Use geometric mean of C* as initial guess
        c_star_mean = np.exp(np.mean(np.log(c_stars + 1e-30)))
        c_oa = self.seed_mass + max(c_star_mean * 0.1, 1.0)

        # Iterative solution
        for iteration in range(max_iterations):
            # Calculate F_p for all species
            f_p = c_oa / (c_oa + c_stars)

            # Calculate new C_OA
            c_particle = c_totals * f_p
            c_oa_new = self.seed_mass + np.sum(c_particle)

            # Check convergence
            if c_oa_new > 0:
                relative_change = abs(c_oa_new - c_oa) / c_oa_new
            else:
                relative_change = abs(c_oa_new - c_oa)

            c_oa = c_oa_new

            if relative_change < tolerance:
                logger.debug(f"Partitioning converged after {iteration + 1} iterations")
                break
        else:
            logger.warning(f"Partitioning did not converge after {max_iterations} iterations")

        # Final calculation with converged C_OA
        f_p = c_oa / (c_oa + c_stars)
        c_particle = c_totals * f_p
        c_gas = c_totals * (1 - f_p)

        # Build results DataFrame
        results = []
        for i, sp in enumerate(species_data):
            results.append({
                'code': sp['code'],
                'formula': sp.get('formula', ''),
                'mw': sp['mw'],
                'log_c_star': np.log10(c_stars[i] + 1e-30),
                'c_star': c_stars[i],
                'c_total': c_totals[i],
                'f_p': f_p[i],
                'c_particle': c_particle[i],
                'c_gas': c_gas[i],
                'c_oa_equilibrium': c_oa
            })

        return pd.DataFrame(results)


# ==============================================================================
# Time Series Processor
# ==============================================================================

# Physical constants for unit conversion
AVOGADRO = 6.022e23  # molecules/mol


def _process_netcdf_concentrations(nc_path: str, species_dict: Dict[str, SpeciesData],
                                   pvap_data: Dict[str, float],
                                   temperature_k: float = 298.0,
                                   seed_mass: float = 0.0) -> Optional[pd.DataFrame]:
    """
    Process concentration data from netCDF file.

    The netCDF file contains:
    - Species: array of species names (with G prefix for gas-phase)
    - time: time array in seconds
    - concentration: 2D array [time, species] in molecules/cm³

    Returns a DataFrame with partitioning calculations for the final timestep.
    """
    try:
        import netCDF4
    except ImportError:
        logger.warning("netCDF4 not available - cannot process netCDF files")
        return None

    try:
        with netCDF4.Dataset(nc_path, 'r') as nc:
            # Get species names
            species_raw = nc.variables['Species'][:]
            species_names = [''.join(s.astype(str)).strip() for s in species_raw]

            # Get concentration data
            conc_data = nc.variables['concentration'][:]
            time_data = nc.variables['time'][:]

            # Use final timestep for partitioning (most relevant for SOA)
            final_conc = conc_data[-1, :]
            final_time = time_data[-1]

            logger.info(f"Processing netCDF: {len(species_names)} species, final time = {final_time}s")

            calculator = PartitioningCalculator(temperature_k, seed_mass)

            # Convert concentrations to µg/m³ and calculate partitioning
            species_data = []

            for i, sp_name in enumerate(species_names):
                conc_molec_cm3 = final_conc[i]
                if conc_molec_cm3 <= 0:
                    continue

                # Strip G prefix for lookup
                lookup_name = sp_name[1:] if sp_name.startswith('G') else sp_name
                sp_info = species_dict.get(lookup_name) or species_dict.get(sp_name)

                if sp_info is None:
                    continue

                # Skip inorganics for partitioning (they don't partition to aerosol)
                if sp_info.n_carbon == 0 and lookup_name not in ['CO', 'CO2', 'CH2O', 'CH3OOH', 'CH3O2']:
                    continue

                # Convert molecules/cm³ to µg/m³
                # µg/m³ = (molec/cm³) * (mol/molec) * (g/mol) * (10^6 cm³/m³) * (10^6 µg/g)
                # µg/m³ = (molec/cm³) * MW / AVOGADRO * 10^12
                mw = sp_info.mw if sp_info.mw > 0 else 100.0  # Default MW if not available
                conc_ug_m3 = conc_molec_cm3 * mw / AVOGADRO * 1e12

                # Get vapor pressure
                log_pvap = pvap_data.get(lookup_name) or pvap_data.get(sp_name)
                if log_pvap is None:
                    log_pvap = VaporPressureReader.estimate_from_mw(
                        mw, sp_info.n_oxygen, sp_info.n_nitrogen
                    )

                c_star = calculator.calculate_c_star(mw, log_pvap)

                species_data.append({
                    'code': lookup_name,
                    'formula': sp_info.formula,
                    'mw': mw,
                    'n_carbon': sp_info.n_carbon,
                    'n_oxygen': sp_info.n_oxygen,
                    'n_nitrogen': sp_info.n_nitrogen,
                    'log_pvap': log_pvap,
                    'c_star': c_star,
                    'concentration': conc_ug_m3
                })

            if not species_data:
                logger.warning("No valid species data extracted from netCDF")
                return None

            logger.info(f"Extracted {len(species_data)} organic species for partitioning")

            # Calculate partitioning using the dynamic algorithm
            df_result = calculator.calculate_partitioning_dynamic(species_data)
            df_result['time'] = final_time

            # Add carbon and nitrogen counts
            code_to_info = {sd['code']: sd for sd in species_data}
            df_result['n_carbon'] = df_result['code'].map(lambda x: code_to_info.get(x, {}).get('n_carbon', 0))
            df_result['n_oxygen'] = df_result['code'].map(lambda x: code_to_info.get(x, {}).get('n_oxygen', 0))
            df_result['n_nitrogen'] = df_result['code'].map(lambda x: code_to_info.get(x, {}).get('n_nitrogen', 0))

            return df_result

    except Exception as e:
        logger.error(f"Error processing netCDF file {nc_path}: {e}")
        import traceback
        traceback.print_exc()
        return None


def process_time_series(output_dir: str, species_dict: Dict[str, SpeciesData],
                        pvap_data: Dict[str, float],
                        temperature_k: float = 298.0,
                        seed_mass: float = 0.0) -> Optional[pd.DataFrame]:
    """
    Process time series data from box model output.

    Reads concentration data and calculates partitioning at each time step
    using the dynamic algorithm.
    """
    # First try netCDF files (primary source of concentration data)
    nc_files = glob.glob(os.path.join(output_dir, "*.nc"))
    if nc_files:
        result = _process_netcdf_concentrations(nc_files[0], species_dict, pvap_data, temperature_k, seed_mass)
        if result is not None:
            return result

    # Fallback to legacy concentration files
    conc_files = glob.glob(os.path.join(output_dir, "*.conc"))
    conc_files.extend(glob.glob(os.path.join(output_dir, "*.out")))
    conc_files.extend(glob.glob(os.path.join(output_dir, "concentration*.dat")))

    if not conc_files:
        logger.warning("No concentration files found for time series processing")
        return None

    calculator = PartitioningCalculator(temperature_k, seed_mass)

    all_results = []

    for conc_file in conc_files:
        try:
            df = _read_concentration_file(conc_file)
            if df is not None and len(df) > 0:
                for idx, row in df.iterrows():
                    time = row.get('time', idx)

                    # Prepare species data for this timestep
                    species_data = []
                    for col in df.columns:
                        if col == 'time':
                            continue

                        conc = row[col]
                        if conc <= 0:
                            continue

                        # Look up species info
                        sp_info = species_dict.get(col) or species_dict.get(f'G{col}')
                        if sp_info is None:
                            continue

                        # Get vapor pressure
                        log_pvap = pvap_data.get(col) or pvap_data.get(sp_info.code)
                        if log_pvap is None:
                            log_pvap = VaporPressureReader.estimate_from_mw(
                                sp_info.mw, sp_info.n_oxygen, sp_info.n_nitrogen
                            )

                        c_star = calculator.calculate_c_star(sp_info.mw, log_pvap)

                        species_data.append({
                            'code': col,
                            'formula': sp_info.formula,
                            'mw': sp_info.mw,
                            'c_star': c_star,
                            'concentration': conc
                        })

                    if species_data:
                        result = calculator.calculate_partitioning_dynamic(species_data)
                        result['time'] = time
                        all_results.append(result)

        except Exception as e:
            logger.warning(f"Error processing {conc_file}: {e}")
            continue

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return None


def _read_concentration_file(filepath: str) -> Optional[pd.DataFrame]:
    """Read a concentration output file."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Find header line
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.strip() and not line.strip().startswith('!'):
                parts = FortranOutputParser.smart_split(line)
                if len(parts) > 2 and any(p.isalpha() for p in parts):
                    # This looks like a header
                    header_line = parts
                    data_start = i + 1
                    break

        if header_line is None:
            return None

        # Read data
        data = []
        for line in lines[data_start:]:
            parts = FortranOutputParser.smart_split(line.strip())
            if len(parts) >= len(header_line):
                try:
                    row = [float(p) for p in parts[:len(header_line)]]
                    data.append(row)
                except ValueError:
                    continue

        if data:
            return pd.DataFrame(data, columns=header_line)
        return None

    except Exception as e:
        logger.debug(f"Could not read {filepath}: {e}")
        return None


def _add_data_quality_note(fig: plt.Figure, note: Optional[str]) -> None:
    """Add a data-quality annotation to plots when synthetic assumptions are used."""
    if not note:
        return
    fig.text(
        0.5, 0.01, note,
        ha='center', va='bottom', fontsize=9, color='darkred'
    )


def _add_synthetic_watermark(fig: plt.Figure, is_synthetic: bool = False,
                             is_surrogate: bool = False,
                             pvap_quality: str = 'good') -> None:
    """
    Add prominent watermark to plots when data quality is compromised.

    This ensures users cannot mistake synthetic/surrogate/degraded data for
    real simulation results.

    Args:
        fig: Matplotlib figure
        is_synthetic: True if concentrations are synthetic (not from simulation)
        is_surrogate: True if data is from baseline surrogate
        pvap_quality: 'good', 'moderate', 'degraded', or 'critical'
    """
    watermark_text = None
    watermark_color = 'lightgray'
    alpha = 0.3

    if is_synthetic:
        watermark_text = "SYNTHETIC DATA"
        watermark_color = 'red'
        alpha = 0.15
    elif is_surrogate:
        watermark_text = "SURROGATE DATA"
        watermark_color = 'orange'
        alpha = 0.15
    elif pvap_quality == 'critical':
        watermark_text = "ESTIMATED PVAP"
        watermark_color = 'red'
        alpha = 0.12
    elif pvap_quality == 'degraded':
        watermark_text = "REDUCED ACCURACY"
        watermark_color = 'darkorange'
        alpha = 0.10

    if watermark_text:
        # Add diagonal watermark across the plot
        fig.text(
            0.5, 0.5, watermark_text,
            fontsize=40, color=watermark_color,
            ha='center', va='center',
            rotation=30, alpha=alpha,
            transform=fig.transFigure,
            fontweight='bold',
            zorder=1000
        )


def _get_plot_metadata(data_quality: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract plot metadata from data quality dict for watermarking decisions.

    Returns:
        Dict with is_synthetic, is_surrogate, pvap_quality keys
    """
    return {
        'is_synthetic': data_quality.get('is_synthetic', False),
        'is_surrogate': data_quality.get('source', '').startswith('baseline') or
                        data_quality.get('is_surrogate', False),
        'pvap_quality': data_quality.get('pvap_quality', 'good')
    }


def _validate_plot_df(
    df: pd.DataFrame,
    required_cols: List[str],
    plot_name: str,
    results: Dict[str, Any]
) -> tuple[bool, Optional[str]]:
    """Validate DataFrame for plot generation and record warnings."""
    if df is None or len(df) == 0:
        reason = "empty dataset"
        results['warnings'].append(f"Skipping {plot_name}: {reason}")
        return False, reason

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        reason = f"missing columns {missing}"
        results['warnings'].append(f"Skipping {plot_name}: {reason}")
        return False, reason

    # Check for any finite data in required columns
    for col in required_cols:
        series = pd.to_numeric(df[col], errors='coerce')
        if series.dropna().empty:
            reason = f"invalid column {col}"
            results['warnings'].append(f"Skipping {plot_name}: {reason}")
            return False, reason

    return True, None


# ==============================================================================
# Plotting Functions
# ==============================================================================

def generate_volatility_distribution_plot(
    df: pd.DataFrame,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate volatility distribution (VBS) plot."""
    if df is None or len(df) == 0:
        return

    # Ensure c_total exists
    df = df.copy()
    if 'c_total' not in df.columns:
        df['c_total'] = df.get('c_particle', 0) + df.get('c_gas', 0)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bin by log10(C*)
    log_c_star = df['log_c_star'].clip(-5, 10)

    # Create bins
    bins = np.arange(-5, 11, 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Sum mass in each bin
    df['bin'] = pd.cut(log_c_star, bins=bins, labels=bin_centers)
    binned = df.groupby('bin', observed=False).agg({
        'c_particle': 'sum',
        'c_gas': 'sum',
        'c_total': 'sum'
    }).reset_index()

    # Plot stacked bars
    width = 0.7
    ax.bar(binned['bin'].astype(float), binned['c_particle'],
           width=width, label='Particle Phase', color='steelblue', alpha=0.8)
    ax.bar(binned['bin'].astype(float), binned['c_gas'],
           width=width, bottom=binned['c_particle'],
           label='Gas Phase', color='lightcoral', alpha=0.8)

    ax.set_xlabel('log$_{10}$(C*) [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_ylabel('Mass Concentration [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_title('Volatility Basis Set Distribution', fontsize=14)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_soa_yield_plot(
    df: pd.DataFrame,
    output_path: str,
    voc_name: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate SOA yield curve plot."""
    if df is None or len(df) == 0:
        return

    # Ensure c_total exists
    df = df.copy()
    if 'c_total' not in df.columns:
        df['c_total'] = df.get('c_particle', 0) + df.get('c_gas', 0)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Group by time if available
    if 'time' in df.columns:
        time_groups = df.groupby('time').agg({
            'c_particle': 'sum',
            'c_total': 'sum',
            'c_oa_equilibrium': 'first'
        }).reset_index()

        ax.scatter(time_groups['c_oa_equilibrium'],
                  time_groups['c_particle'] / (time_groups['c_total'] + 1e-10),
                  c=time_groups['time'], cmap='viridis', s=50, alpha=0.7)

        cbar = plt.colorbar(ax.collections[0])
        cbar.set_label('Time', fontsize=10)
    else:
        # Plot all points without time coloring
        c_oa = df['c_oa_equilibrium']
        soa_yield = df['c_particle'] / (df['c_total'] + 1e-10)
        ax.scatter(c_oa, soa_yield, s=50, c='steelblue', alpha=0.6, label='Data Points')

    ax.set_xlabel('C$_{OA}$ [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_ylabel('SOA Yield (M$_{OA}$ / $\\Delta$VOC)', fontsize=12)
    ax.set_title(f'SOA Yield Curve: {voc_name}', fontsize=14)
    ax.grid(alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_partitioning_summary_plot(
    df: pd.DataFrame,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate partitioning factor summary plot."""
    if df is None or len(df) == 0:
        return

    # Ensure c_total exists
    df = df.copy()
    if 'c_total' not in df.columns:
        df['c_total'] = df.get('c_particle', 0) + df.get('c_gas', 0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: F_p vs log(C*)
    ax1 = axes[0]
    c_total_positive = df['c_total'].clip(lower=1e-10)
    scatter = ax1.scatter(df['log_c_star'], df['f_p'],
                          c=c_total_positive, cmap='plasma',
                          s=30, alpha=0.6, norm=matplotlib.colors.LogNorm())

    ax1.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='F$_p$ = 0.5')
    ax1.set_xlabel('log$_{10}$(C*) [$\\mu$g/m$^3$]', fontsize=12)
    ax1.set_ylabel('Particle Fraction (F$_p$)', fontsize=12)
    ax1.set_title('Partitioning by Volatility', fontsize=12)
    ax1.set_ylim(0, 1)
    ax1.legend()

    cbar1 = plt.colorbar(scatter, ax=ax1)
    cbar1.set_label('Total Concentration', fontsize=10)

    # Right: Pie chart of phase distribution
    ax2 = axes[1]
    total_particle = df['c_particle'].sum()
    total_gas = df['c_gas'].sum()

    if total_particle + total_gas > 0:
        sizes = [total_particle, total_gas]
        labels = [f'Particle\n{total_particle:.2f} $\\mu$g/m$^3$',
                 f'Gas\n{total_gas:.2f} $\\mu$g/m$^3$']
        colors = ['steelblue', 'lightcoral']
        ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                startangle=90, explode=(0.05, 0))
        ax2.set_title('Phase Distribution', fontsize=12)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_van_krevelen_plot(
    df: pd.DataFrame,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate Van Krevelen diagram (H:C vs O:C)."""
    if df is None or len(df) == 0:
        return

    # Calculate H:C and O:C ratios
    if 'n_carbon' not in df.columns or 'n_oxygen' not in df.columns:
        return

    df_plot = df[(df['n_carbon'] > 0)].copy()
    if len(df_plot) == 0:
        return

    # Estimate H from formula (MW - contributions from C, O, N)
    # Using integer atomic weights consistent with GECKO mechanism definitions
    n_nitrogen = df_plot['n_nitrogen'] if 'n_nitrogen' in df_plot.columns else 0
    df_plot['n_hydrogen'] = (df_plot['mw'] - 12.0*df_plot['n_carbon'] - 16.0*df_plot['n_oxygen'] - 14.0*n_nitrogen) / 1.008
    df_plot['n_hydrogen'] = df_plot['n_hydrogen'].clip(0, None)

    df_plot['hc_ratio'] = df_plot['n_hydrogen'] / df_plot['n_carbon']
    df_plot['oc_ratio'] = df_plot['n_oxygen'] / df_plot['n_carbon']

    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by particle mass
    if 'c_particle' in df_plot.columns and df_plot['c_particle'].sum() > 0:
        scatter = ax.scatter(df_plot['oc_ratio'], df_plot['hc_ratio'],
                            c=df_plot['c_particle'], cmap='viridis',
                            s=30, alpha=0.6, norm=matplotlib.colors.LogNorm(vmin=1e-10))
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Particle Mass [$\\mu$g/m$^3$]', fontsize=10)
    else:
        ax.scatter(df_plot['oc_ratio'], df_plot['hc_ratio'],
                  c='steelblue', s=30, alpha=0.6)

    # Add reference lines for oxidation pathways
    oc = np.linspace(0, 2, 100)
    ax.plot(oc, 2 - oc, 'k--', alpha=0.3, label='Carboxylic acid addition')
    ax.plot(oc, 2 - 2*oc, 'k:', alpha=0.3, label='Ketone/aldehyde formation')

    ax.set_xlabel('O:C Ratio', fontsize=12)
    ax.set_ylabel('H:C Ratio', fontsize=12)
    ax.set_title('Van Krevelen Diagram', fontsize=14)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2.5)
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_top10_species_plot(
    df: pd.DataFrame,
    output_path: str,
    phase: str = 'gas',
    data_quality: Optional[Dict[str, Any]] = None,
    data_quality_note: Optional[str] = None
):
    """Generate top 10 species bar chart for gas or particle phase."""
    if df is None or len(df) == 0:
        return

    col = 'c_gas' if phase == 'gas' else 'c_particle'
    if col not in df.columns:
        return

    # Filter out inorganics and non-organic species
    INORGANIC_CODES = {
        'H2', 'H2O', 'H2O2', 'HNO2', 'HNO3', 'HNO4', 'HO', 'HO2',
        'N2O5', 'NO', 'NO2', 'NO3', 'O3P', 'O1D', 'O2', 'O3',
        'SO2', 'SULF', 'XCLOST', 'CO', 'CO2', 'N2', 'CH4', 'CH2O', 'CH3OOH',
        'GH2', 'GH2O', 'GH2O2', 'GHNO2', 'GHNO3', 'GHNO4', 'GHO', 'GHO2',
        'GN2O5', 'GNO', 'GNO2', 'GNO3', 'GO3P', 'GO1D', 'GO2', 'GO3',
        'GSO2', 'GSULF', 'GXCLOST', 'GCO', 'GCO2', 'GCH4', 'GCH2O', 'GCH3OOH'
    }

    df_filtered = df.copy()

    # Filter by code (exclude inorganics and small molecules)
    if 'code' in df_filtered.columns:
        df_filtered = df_filtered[~df_filtered['code'].isin(INORGANIC_CODES)]

    # Filter by n_carbon (must have more than 1 carbon for meaningful organic species)
    if 'n_carbon' in df_filtered.columns:
        df_filtered = df_filtered[df_filtered['n_carbon'] > 1]

    # Filter out zero/very small values
    df_filtered = df_filtered[df_filtered[col] > 1e-20]

    if len(df_filtered) == 0:
        logger.warning(f"No valid {phase} species after filtering")
        return

    # Get top 10 by concentration
    df_sorted = df_filtered.nlargest(10, col)

    if len(df_sorted) == 0 or df_sorted[col].sum() == 0:
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    bars = ax.barh(range(len(df_sorted)), df_sorted[col].values,
                   color='steelblue' if phase == 'gas' else 'coral', alpha=0.8)

    # Labels: use code if formula is too long
    labels = []
    for _, row in df_sorted.iterrows():
        label = row.get('code', row.get('formula', 'Unknown'))
        if len(label) > 20:
            label = label[:17] + '...'
        labels.append(label)

    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    ax.set_xlabel('Concentration [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_title(f'Top 10 {"Gas" if phase == "gas" else "Particle"}-Phase Species', fontsize=14)
    ax.grid(axis='x', alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_trajectory_plot(
    df: pd.DataFrame,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate NC vs O:C trajectory plot."""
    if df is None or len(df) == 0:
        return

    if 'n_carbon' not in df.columns or 'n_oxygen' not in df.columns:
        return

    df_plot = df[(df['n_carbon'] > 0)].copy()
    if len(df_plot) == 0:
        return

    df_plot['oc_ratio'] = df_plot['n_oxygen'] / df_plot['n_carbon']

    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by particle fraction
    if 'f_p' in df_plot.columns:
        scatter = ax.scatter(df_plot['n_carbon'], df_plot['oc_ratio'],
                            c=df_plot['f_p'], cmap='RdYlBu_r',
                            s=30, alpha=0.6, vmin=0, vmax=1)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Particle Fraction (F$_p$)', fontsize=10)
    else:
        ax.scatter(df_plot['n_carbon'], df_plot['oc_ratio'],
                  c='steelblue', s=30, alpha=0.6)

    ax.set_xlabel('Carbon Number (n$_C$)', fontsize=12)
    ax.set_ylabel('O:C Ratio', fontsize=12)
    ax.set_title('Carbon Number vs. Oxidation Trajectory', fontsize=14)
    ax.grid(alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_soa_mass_plot(
    df: pd.DataFrame,
    output_path: str,
    voc_name: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate SOA mass time series plot."""
    if df is None or len(df) == 0:
        return

    if 'time' not in df.columns:
        # Single time point - just show total
        fig, ax = plt.subplots(figsize=(10, 6))
        total_soa = df['c_particle'].sum() if 'c_particle' in df.columns else 0
        ax.bar([voc_name], [total_soa], color='steelblue')
        ax.set_ylabel('SOA Mass [$\\mu$g/m$^3$]', fontsize=12)
        ax.set_title(f'Total SOA Mass: {voc_name}', fontsize=14)
        # Add watermark if data quality is compromised
        if data_quality:
            meta = _get_plot_metadata(data_quality)
            _add_synthetic_watermark(fig, **meta)
        _add_data_quality_note(fig, data_quality_note)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()
        return

    # Time series
    time_groups = df.groupby('time').agg({
        'c_particle': 'sum',
        'c_gas': 'sum'
    }).reset_index()

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.fill_between(time_groups['time'], 0, time_groups['c_particle'],
                    alpha=0.7, color='steelblue', label='Particle Phase')
    ax.plot(time_groups['time'], time_groups['c_particle'],
            color='darkblue', linewidth=2)

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('SOA Mass [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_title(f'SOA Mass Evolution: {voc_name}', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_functional_group_plot(
    df: pd.DataFrame,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate functional group distribution plot."""
    if df is None or len(df) == 0:
        return

    # Define functional groups based on composition
    groups = {
        'Alcohols (-OH)': 0,
        'Carbonyls (C=O)': 0,
        'Carboxylic acids (-COOH)': 0,
        'Nitrates (-ONO2)': 0,
        'Peroxides (-OOH)': 0,
        'Ethers (-O-)': 0,
        'Other': 0
    }

    # Count based on oxygen and nitrogen content
    for _, row in df.iterrows():
        n_o = row.get('n_oxygen', 0)
        n_n = row.get('n_nitrogen', 0)
        mass = row.get('c_particle', 0) + row.get('c_gas', 0)

        if n_n > 0:
            groups['Nitrates (-ONO2)'] += mass * n_n
        if n_o > n_n * 3:  # Remaining oxygen
            remaining_o = n_o - n_n * 3
            # Simple heuristic distribution
            groups['Alcohols (-OH)'] += mass * remaining_o * 0.3
            groups['Carbonyls (C=O)'] += mass * remaining_o * 0.3
            groups['Peroxides (-OOH)'] += mass * remaining_o * 0.2
            groups['Ethers (-O-)'] += mass * remaining_o * 0.2

    # Filter out zeros
    groups = {k: v for k, v in groups.items() if v > 0}

    if not groups:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
    bars = ax.bar(range(len(groups)), list(groups.values()), color=colors)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(list(groups.keys()), rotation=45, ha='right')
    ax.set_ylabel('Estimated Mass [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_title('Functional Group Distribution (Estimated)', fontsize=14)
    ax.grid(axis='y', alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_osc_volatility_plot(
    df: pd.DataFrame,
    output_path: str,
    label: str = '',
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate oxidation state of carbon (OSc) vs volatility plot."""
    if df is None or len(df) == 0:
        return

    if 'n_carbon' not in df.columns or 'n_oxygen' not in df.columns:
        return

    df_plot = df[(df['n_carbon'] > 0)].copy()
    if len(df_plot) == 0:
        return

    # Calculate OSc = 2*O:C - H:C (approximation)
    # For organics: OSc ≈ 2*(O/C) - (H/C)
    df_plot['oc_ratio'] = df_plot['n_oxygen'] / df_plot['n_carbon']
    n_nitrogen = df_plot['n_nitrogen'] if 'n_nitrogen' in df_plot.columns else 0
    df_plot['n_hydrogen'] = (df_plot['mw'] - 12*df_plot['n_carbon'] - 16*df_plot['n_oxygen'] - 14*n_nitrogen) / 1.008
    df_plot['n_hydrogen'] = df_plot['n_hydrogen'].clip(0, None)
    df_plot['hc_ratio'] = df_plot['n_hydrogen'] / df_plot['n_carbon']
    df_plot['osc'] = 2 * df_plot['oc_ratio'] - df_plot['hc_ratio']

    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by total mass
    if 'c_particle' in df_plot.columns:
        total_mass = df_plot['c_particle'] + df_plot.get('c_gas', 0)
        scatter = ax.scatter(df_plot['log_c_star'], df_plot['osc'],
                            c=total_mass, cmap='plasma',
                            s=30, alpha=0.6, norm=matplotlib.colors.LogNorm(vmin=1e-10))
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Total Mass [$\\mu$g/m$^3$]', fontsize=10)
    else:
        ax.scatter(df_plot['log_c_star'], df_plot['osc'],
                  c='steelblue', s=30, alpha=0.6)

    # Add reference regions
    ax.axhspan(-2, 0, alpha=0.1, color='green', label='Fresh SOA')
    ax.axhspan(0, 2, alpha=0.1, color='orange', label='Aged SOA')

    ax.set_xlabel('log$_{10}$(C*) [$\\mu$g/m$^3$]', fontsize=12)
    ax.set_ylabel('OS$_C$ (Oxidation State of Carbon)', fontsize=12)
    ax.set_title(f'OSc vs Volatility {label}', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)

    # Add watermark if data quality is compromised
    if data_quality:
        meta = _get_plot_metadata(data_quality)
        _add_synthetic_watermark(fig, **meta)

    _add_data_quality_note(fig, data_quality_note)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def generate_main_simulation_plot(
    output_dir: str,
    output_path: str,
    data_quality_note: Optional[str] = None,
    data_quality: Optional[Dict[str, Any]] = None
):
    """Generate main simulation overview plot from netCDF or time series data."""
    import os

    # Try to read netCDF file
    nc_files = [f for f in os.listdir(output_dir) if f.endswith('.nc')]
    if nc_files:
        try:
            import netCDF4
            nc_path = os.path.join(output_dir, nc_files[0])
            with netCDF4.Dataset(nc_path) as nc:
                # Extract time variable
                time = nc.variables.get('time', nc.variables.get('TIME'))
                if time is not None:
                    time_data = time[:]

                    # Get species names from Species variable
                    species_names = []
                    if 'Species' in nc.variables:
                        species_raw = nc.variables['Species'][:]
                        species_names = [''.join(s.astype(str)).strip() for s in species_raw]

                    # Get concentration data
                    conc_data = None
                    if 'concentration' in nc.variables:
                        conc_data = nc.variables['concentration'][:]

                    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

                    # Plot key oxidants - try both with and without G prefix
                    oxidant_map = {
                        'O3': ['GO3', 'O3', 'ozone'],
                        'OH': ['GHO', 'OH', 'HO'],
                        'NO': ['GNO', 'NO'],
                        'NO2': ['GNO2', 'NO2']
                    }

                    for i, (display_name, var_names) in enumerate(oxidant_map.items()):
                        ax = axes[i // 2, i % 2]
                        plotted = False

                        # First try direct variable access
                        for vn in var_names:
                            if vn in nc.variables:
                                ax.plot(time_data, nc.variables[vn][:], linewidth=1.5, color='steelblue')
                                ax.set_xlabel('Time [s]')
                                ax.set_ylabel(f'{display_name} [molec/cm³]')
                                ax.set_title(display_name)
                                ax.grid(alpha=0.3)
                                plotted = True
                                break

                        # If not found, try species array
                        if not plotted and species_names and conc_data is not None:
                            for vn in var_names:
                                if vn in species_names:
                                    idx = species_names.index(vn)
                                    ax.plot(time_data, conc_data[:, idx], linewidth=1.5, color='steelblue')
                                    ax.set_xlabel('Time [s]')
                                    ax.set_ylabel(f'{display_name} [molec/cm³]')
                                    ax.set_title(display_name)
                                    ax.grid(alpha=0.3)
                                    plotted = True
                                    break

                        if not plotted:
                            ax.text(0.5, 0.5, f'{display_name}\nnot available',
                                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
                            ax.set_frame_on(False)
                            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

                    plt.suptitle('Key Oxidant Time Series', fontsize=14, fontweight='bold')
                    # Add watermark if data quality is compromised
                    if data_quality:
                        meta = _get_plot_metadata(data_quality)
                        _add_synthetic_watermark(fig, **meta)
                    _add_data_quality_note(fig, data_quality_note)
                    plt.tight_layout()
                    plt.savefig(output_path, dpi=300)
                    plt.close()
                    return
        except Exception as e:
            logger.warning(f"Could not read netCDF: {e}")

    # Fallback: create a placeholder plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, 'Time series data not available\nRun box model for detailed simulation results',
           ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # Add watermark for placeholder
    _add_synthetic_watermark(fig, is_synthetic=True)
    _add_data_quality_note(fig, data_quality_note)
    plt.savefig(output_path, dpi=300)
    plt.close()


# ==============================================================================
# Main Processing Function
# ==============================================================================

def run_postprocessing(output_dir: str, voc_name: str,
                       temperature_k: float = 298.0,
                       seed_mass: float = 0.0,
                       allow_synthetic: bool = True) -> Dict[str, Any]:
    """
    Run complete post-processing pipeline.

    Args:
        output_dir: Directory containing GECKO-A output files
        voc_name: Name of the VOC being processed
        temperature_k: Temperature for calculations
        seed_mass: Initial seed aerosol mass (ug/m3)

    Returns:
        Dict with processing results and output paths
    """
    results = {
        'status': 'success',
        'warnings': [],
        'outputs': {},
        'data_quality': {
            'source': 'simulation',
            'is_synthetic': False,
            'assumptions': [],
            'data_points': 0,
            'time_points': 0
        }
    }

    logger.info(f"Running post-processing for {voc_name} in {output_dir}")

    try:
        # 1. Parse dictionary
        dict_path = os.path.join(output_dir, "dictionary.out")
        species_dict = parse_dictionary_robust(dict_path)

        if not species_dict:
            results['warnings'].append("No species data found in dictionary")
            return results

        # 2. Read vapor pressures
        pvap_data = VaporPressureReader.read_pvap_file(output_dir)
        if not pvap_data:
            # CRITICAL WARNING: MW-based Pvap estimation can be off by 2+ orders of magnitude
            results['warnings'].append(
                "⚠️ CRITICAL: Vapor pressure data file (pvap.nan.dat) not found. "
                "Using MW-based estimates for ALL species. Partitioning accuracy may be "
                "reduced by 2+ orders of magnitude. Run GECKO-A generator to produce "
                "calculated vapor pressures using Nannoolal SAR."
            )
            results['data_quality']['assumptions'].append(
                "ALL Pvap values estimated using MW correlation (pvap file missing) - "
                "partitioning accuracy significantly reduced"
            )
            results['data_quality']['pvap_source'] = 'mw_estimation_only'
            results['data_quality']['pvap_accuracy_warning'] = True

        # 3. Calculate C* for all species
        calculator = PartitioningCalculator(temperature_k, seed_mass)

        species_data = []
        pvap_estimated_count = 0
        pvap_total_count = 0
        for code, sp in species_dict.items():
            if code.startswith('G') and code[1:] in species_dict:
                continue  # Skip duplicates

            log_pvap = pvap_data.get(code)
            if log_pvap is None:
                log_pvap = VaporPressureReader.estimate_from_mw(
                    sp.mw, sp.n_oxygen, sp.n_nitrogen
                )
                if pvap_data:  # Only warn if we have some pvap data
                    results['warnings'].append(f"Using estimated Pvap for {code}")
                pvap_estimated_count += 1
            pvap_total_count += 1

            c_star = calculator.calculate_c_star(sp.mw, log_pvap)

            species_data.append({
                'code': code,
                'formula': sp.formula,
                'mw': sp.mw,
                'n_carbon': sp.n_carbon,
                'n_oxygen': sp.n_oxygen,
                'log_pvap': log_pvap,
                'c_star': c_star,
                'concentration': 1.0  # Default; will be overwritten if conc data exists
            })

        # 4. Create base aerosol data (from dictionary only)
        df_base = pd.DataFrame(species_data)
        df_base['log_c_star'] = np.log10(df_base['c_star'] + 1e-30)
        c_oa_assumed = max(seed_mass, 10.0)
        df_base['f_p'] = c_oa_assumed / (c_oa_assumed + df_base['c_star'])  # Static estimate
        df_base['c_particle'] = df_base['concentration'] * df_base['f_p']
        df_base['c_gas'] = df_base['concentration'] * (1 - df_base['f_p'])
        df_base['c_oa_equilibrium'] = c_oa_assumed

        # 5. Try to process time series if available
        df_timeseries = process_time_series(
            output_dir, species_dict, pvap_data, temperature_k, seed_mass
        )

        data_quality_note = None
        # Use time series data if available, otherwise base data
        if df_timeseries is not None and len(df_timeseries) > 0:
            df_final = df_timeseries
        else:
            if not allow_synthetic:
                results['status'] = 'error'
                results['error'] = 'No simulation concentration data found; synthetic fallback disabled.'
                results['warnings'].append(results['error'])
                return results

            results['data_quality']['source'] = 'synthetic_dictionary'
            results['data_quality']['is_synthetic'] = True
            results['data_quality']['assumptions'] = [
                'Concentrations set to 1.0 µg/m³ for all species',
                f'Static partitioning with C_OA = {c_oa_assumed:.2f} µg/m³'
            ]
            results['warnings'].append(
                'No simulation concentration data found. Using synthetic base dataset from dictionary with assumed concentrations.'
            )
            data_quality_note = (
                f'Synthetic data: concentrations assumed; C_OA={c_oa_assumed:.2f} µg/m³. '
                'Run Box Model for real results.'
            )
            df_final = df_base

        results['data_quality']['data_points'] = int(len(df_final))
        if 'time' in df_final.columns:
            results['data_quality']['time_points'] = int(df_final['time'].nunique())

        if pvap_total_count > 0:
            results['data_quality']['pvap_estimated_count'] = int(pvap_estimated_count)
            results['data_quality']['pvap_total_count'] = int(pvap_total_count)
            pvap_estimated_fraction = pvap_estimated_count / pvap_total_count

            # Track vapor pressure data quality with severity levels
            results['data_quality']['pvap_estimated_fraction'] = float(pvap_estimated_fraction)

            if pvap_estimated_fraction >= 0.9:
                # CRITICAL: Almost all Pvap estimated
                results['warnings'].append(
                    f"⚠️ CRITICAL: {pvap_estimated_fraction*100:.0f}% of species "
                    f"({pvap_estimated_count}/{pvap_total_count}) used MW-based Pvap estimates. "
                    "Partitioning results may be inaccurate by 2+ orders of magnitude."
                )
                results['data_quality']['assumptions'].append(
                    f"CRITICAL: {pvap_estimated_fraction*100:.0f}% Pvap values estimated - "
                    "partitioning accuracy severely compromised"
                )
                results['data_quality']['pvap_quality'] = 'critical'
            elif pvap_estimated_fraction >= 0.5:
                # HIGH: More than half estimated
                results['warnings'].append(
                    f"⚠️ WARNING: {pvap_estimated_fraction*100:.0f}% of species "
                    f"({pvap_estimated_count}/{pvap_total_count}) used MW-based Pvap estimates. "
                    "Partitioning accuracy reduced."
                )
                results['data_quality']['assumptions'].append(
                    f"High fraction ({pvap_estimated_fraction*100:.0f}%) of Pvap values estimated"
                )
                results['data_quality']['pvap_quality'] = 'degraded'
            elif pvap_estimated_fraction >= 0.2:
                # MODERATE: Significant fraction estimated
                results['warnings'].append(
                    f"Note: {pvap_estimated_fraction*100:.0f}% of species used estimated Pvap"
                )
                results['data_quality']['pvap_quality'] = 'moderate'
            else:
                # GOOD: Most Pvap from calculated values
                results['data_quality']['pvap_quality'] = 'good'

        df_final['data_source'] = results['data_quality']['source']

        # Sanitize invalid values
        df_final = df_final.replace([np.inf, -np.inf], np.nan)
        for col in ['c_particle', 'c_gas', 'c_total']:
            if col in df_final.columns:
                neg_count = (df_final[col] < 0).sum()
                if neg_count > 0:
                    results['warnings'].append(f"Clipped {neg_count} negative values in {col}")
                df_final[col] = df_final[col].clip(lower=0)

        # 6. Save aerosol data
        csv_path = os.path.join(output_dir, "aerosol_data.csv")
        df_final.to_csv(csv_path, index=False)
        results['outputs']['aerosol_data'] = csv_path

        # 7. Generate plots (strict validation gates)
        plot_audit = {}

        vbs_plot = os.path.join(output_dir, "volatility_distribution.png")
        ok, reason = _validate_plot_df(df_final, ['log_c_star', 'c_particle', 'c_gas'], 'VBS plot', results)
        plot_audit['vbs_plot'] = {
            'title': 'Volatility Basis Set Distribution',
            'file': vbs_plot,
            'required_columns': ['log_c_star', 'c_particle', 'c_gas'],
            'units': {
                'log_c_star': 'log10(µg/m³)',
                'c_particle': 'µg/m³',
                'c_gas': 'µg/m³'
            },
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_volatility_distribution_plot(df_final, vbs_plot, data_quality_note, results['data_quality'])
            results['outputs']['vbs_plot'] = vbs_plot

        yield_plot = os.path.join(output_dir, "soa_yield.png")
        ok, reason = _validate_plot_df(df_final, ['c_particle', 'c_total'], 'SOA yield plot', results)
        plot_audit['yield_plot'] = {
            'title': 'SOA Yield',
            'file': yield_plot,
            'required_columns': ['c_particle', 'c_total'],
            'units': {
                'c_particle': 'µg/m³',
                'c_total': 'µg/m³'
            },
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_soa_yield_plot(df_final, yield_plot, voc_name, data_quality_note, results['data_quality'])
            results['outputs']['yield_plot'] = yield_plot

        summary_plot = os.path.join(output_dir, "partitioning_summary.png")
        ok, reason = _validate_plot_df(df_final, ['log_c_star', 'f_p', 'c_particle', 'c_gas'], 'Partitioning summary plot', results)
        plot_audit['summary_plot'] = {
            'title': 'Partitioning Summary',
            'file': summary_plot,
            'required_columns': ['log_c_star', 'f_p', 'c_particle', 'c_gas'],
            'units': {
                'log_c_star': 'log10(µg/m³)',
                'f_p': 'fraction',
                'c_particle': 'µg/m³',
                'c_gas': 'µg/m³'
            },
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_partitioning_summary_plot(df_final, summary_plot, data_quality_note, results['data_quality'])
            results['outputs']['summary_plot'] = summary_plot

        # Additional plots (restored from v1.0)
        vk_plot = os.path.join(output_dir, "van_krevelen_HC_OC.png")
        ok, reason = _validate_plot_df(df_final, ['n_carbon', 'n_oxygen', 'mw'], 'Van Krevelen plot', results)
        plot_audit['van_krevelen'] = {
            'title': 'Van Krevelen (H:C vs O:C)',
            'file': vk_plot,
            'required_columns': ['n_carbon', 'n_oxygen', 'mw'],
            'units': {
                'n_carbon': 'count',
                'n_oxygen': 'count',
                'mw': 'g/mol'
            },
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_van_krevelen_plot(df_final, vk_plot, data_quality_note, results['data_quality'])
            results['outputs']['van_krevelen'] = vk_plot

        top10_gas_plot = os.path.join(output_dir, "top10_gas.png")
        ok, reason = _validate_plot_df(df_final, ['c_gas'], 'Top-10 gas plot', results)
        plot_audit['top10_gas'] = {
            'title': 'Top-10 Gas Species',
            'file': top10_gas_plot,
            'required_columns': ['c_gas'],
            'units': {'c_gas': 'µg/m³'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_top10_species_plot(df_final, top10_gas_plot, phase='gas', data_quality=results['data_quality'], data_quality_note=data_quality_note)
            results['outputs']['top10_gas'] = top10_gas_plot

        top10_particle_plot = os.path.join(output_dir, "top10_particle.png")
        ok, reason = _validate_plot_df(df_final, ['c_particle'], 'Top-10 particle plot', results)
        plot_audit['top10_particle'] = {
            'title': 'Top-10 Particle Species',
            'file': top10_particle_plot,
            'required_columns': ['c_particle'],
            'units': {'c_particle': 'µg/m³'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_top10_species_plot(df_final, top10_particle_plot, phase='particle', data_quality=results['data_quality'], data_quality_note=data_quality_note)
            results['outputs']['top10_particle'] = top10_particle_plot

        trajectory_plot = os.path.join(output_dir, "trajectory_NC_OC.png")
        ok, reason = _validate_plot_df(df_final, ['n_carbon', 'n_oxygen'], 'Trajectory plot', results)
        plot_audit['trajectory'] = {
            'title': 'Carbon Number vs O:C Trajectory',
            'file': trajectory_plot,
            'required_columns': ['n_carbon', 'n_oxygen'],
            'units': {'n_carbon': 'count', 'n_oxygen': 'count'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_trajectory_plot(df_final, trajectory_plot, data_quality_note, results['data_quality'])
            results['outputs']['trajectory'] = trajectory_plot

        soa_mass_plot = os.path.join(output_dir, "soa_mass.png")
        ok, reason = _validate_plot_df(df_final, ['c_particle'], 'SOA mass plot', results)
        plot_audit['soa_mass'] = {
            'title': 'SOA Mass',
            'file': soa_mass_plot,
            'required_columns': ['c_particle'],
            'units': {'c_particle': 'µg/m³'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_soa_mass_plot(df_final, soa_mass_plot, voc_name, data_quality_note, results['data_quality'])
            results['outputs']['soa_mass'] = soa_mass_plot

        fg_plot = os.path.join(output_dir, "functional_group_dist_N.png")
        ok, reason = _validate_plot_df(df_final, ['n_oxygen'], 'Functional group plot', results)
        plot_audit['functional_groups'] = {
            'title': 'Functional Group Distribution',
            'file': fg_plot,
            'required_columns': ['n_oxygen'],
            'units': {'n_oxygen': 'count'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_functional_group_plot(df_final, fg_plot, data_quality_note, results['data_quality'])
            results['outputs']['functional_groups'] = fg_plot

        osc_plot = os.path.join(output_dir, "osc_volatility_1x_Lifetime_Approx.png")
        ok, reason = _validate_plot_df(df_final, ['n_carbon', 'n_oxygen', 'mw', 'log_c_star'], 'OSc volatility plot', results)
        plot_audit['osc_volatility'] = {
            'title': 'OSc vs Volatility',
            'file': osc_plot,
            'required_columns': ['n_carbon', 'n_oxygen', 'mw', 'log_c_star'],
            'units': {
                'n_carbon': 'count',
                'n_oxygen': 'count',
                'mw': 'g/mol',
                'log_c_star': 'log10(µg/m³)'
            },
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': ok,
            'skipped_reason': reason
        }
        if ok:
            generate_osc_volatility_plot(df_final, osc_plot, label='1x Lifetime', data_quality_note=data_quality_note, data_quality=results['data_quality'])
            results['outputs']['osc_volatility'] = osc_plot

        main_plot = os.path.join(output_dir, "plot.png")
        plot_audit['main_plot'] = {
            'title': 'Key Oxidant Time Series',
            'file': main_plot,
            'required_columns': ['time', 'concentration (netCDF)'],
            'units': {'time': 's', 'concentration': 'molec/cm³'},
            'data_source': results['data_quality']['source'],
            'assumptions': results['data_quality'].get('assumptions', []),
            'generated': True,
            'skipped_reason': None
        }
        generate_main_simulation_plot(output_dir, main_plot, data_quality_note, results['data_quality'])
        results['outputs']['main_plot'] = main_plot

        # 8. Calculate summary statistics
        total_soa = df_final['c_particle'].sum() if 'c_particle' in df_final.columns else 0.0
        total_gas = df_final['c_gas'].sum() if 'c_gas' in df_final.columns else 0.0
        total_organic = total_soa + total_gas

        soa_yield = total_soa / total_organic if total_organic > 0 else 0

        results['summary'] = {
            'total_species': len(species_dict),
            'total_soa_mass_ug_m3': float(total_soa),
            'total_gas_mass_ug_m3': float(total_gas),
            'soa_yield': float(soa_yield),
            'mean_c_star': float(df_final['c_star'].mean()) if 'c_star' in df_final.columns else 0.0,
            'median_log_c_star': float(df_final['log_c_star'].median()) if 'log_c_star' in df_final.columns else 0.0,
            'c_oa_equilibrium': float(df_final['c_oa_equilibrium'].mean()) if 'c_oa_equilibrium' in df_final.columns else 0.0,
            'data_quality': results['data_quality'],
            'plot_audit': plot_audit
        }

        # Save summary
        import json
        summary_path = os.path.join(output_dir, "postprocessing_summary.json")
        with open(summary_path, 'w') as f:
            json.dump(results['summary'], f, indent=2)
        results['outputs']['summary'] = summary_path

        plot_audit_path = os.path.join(output_dir, "plot_audit.json")
        with open(plot_audit_path, 'w') as f:
            json.dump(plot_audit, f, indent=2)
        results['outputs']['plot_audit'] = plot_audit_path

        # Data Quality Appendix (JSON + PDF)
        appendix = {
            'voc_name': voc_name,
            'data_quality': results['data_quality'],
            'plot_audit': plot_audit,
            'warnings': results.get('warnings', [])
        }
        appendix_json_path = os.path.join(output_dir, "data_quality_appendix.json")
        with open(appendix_json_path, 'w') as f:
            json.dump(appendix, f, indent=2)
        results['outputs']['data_quality_appendix_json'] = appendix_json_path

        try:
            from matplotlib.backends.backend_pdf import PdfPages

            appendix_pdf_path = os.path.join(output_dir, "data_quality_appendix.pdf")
            with PdfPages(appendix_pdf_path) as pdf:
                fig = plt.figure(figsize=(8.27, 11.69))
                fig.text(0.02, 0.98, f"Data Quality Appendix - {voc_name}", fontsize=14, fontweight='bold')

                lines = [
                    f"Source: {results['data_quality'].get('source')}",
                    f"Synthetic: {results['data_quality'].get('is_synthetic')}",
                    f"Data Points: {results['data_quality'].get('data_points')}",
                    f"Time Points: {results['data_quality'].get('time_points')}",
                    f"Pvap Estimated: {results['data_quality'].get('pvap_estimated_count', 0)} / {results['data_quality'].get('pvap_total_count', 0)}",
                    "",
                    "Assumptions:",
                ]
                assumptions = results['data_quality'].get('assumptions', [])
                if not assumptions:
                    assumptions = ["None"]
                for a in assumptions:
                    lines.append(f"- {a}")
                lines.append("")
                lines.append("Warnings:")
                warnings = results.get('warnings', []) or ["None"]
                for w in warnings:
                    lines.append(f"- {w}")

                y = 0.94
                for line in lines:
                    fig.text(0.02, y, line, fontsize=10)
                    y -= 0.02

                pdf.savefig(fig)
                plt.close(fig)

            results['outputs']['data_quality_appendix_pdf'] = appendix_pdf_path
        except Exception as e:
            results['warnings'].append(f"Failed to create data quality PDF: {e}")

        logger.info(f"Post-processing complete. SOA yield: {soa_yield:.3f}")

    except Exception as e:
        logger.error(f"Post-processing failed: {e}")
        results['status'] = 'error'
        results['error'] = str(e)

    return results


# ==============================================================================
# CLI Entry Point
# ==============================================================================

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python postprocessing.py <output_dir> [voc_name] [temperature_K] [seed_mass]")
        sys.exit(1)

    output_dir = sys.argv[1]
    voc_name = sys.argv[2] if len(sys.argv) > 2 else "unknown"
    temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 298.0
    seed = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0

    logging.basicConfig(level=logging.INFO)

    results = run_postprocessing(output_dir, voc_name, temperature, seed)

    print(f"\nStatus: {results['status']}")
    if results.get('summary'):
        print(f"SOA Mass: {results['summary']['total_soa_mass_ug_m3']:.3f} ug/m3")
        print(f"SOA Yield: {results['summary']['soa_yield']:.4f}")
        print(f"C_OA (equilibrium): {results['summary']['c_oa_equilibrium']:.2f} ug/m3")
