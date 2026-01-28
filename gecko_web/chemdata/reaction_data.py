"""
GECKO-A Reaction Rate Constants and Kinetics Database

This module provides experimentally validated rate constants and branching ratios
for atmospheric oxidation reactions used in GECKO-A.

Data Sources:
- IUPAC Task Group on Atmospheric Chemical Kinetic Data Evaluation
  (https://iupac.aeris-data.fr/)
- NASA/JPL Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies
  (JPL Publication 19-5, 2020)
- MCM (Master Chemical Mechanism) v3.3.1
- Atkinson et al. (2006): Evaluated kinetic and photochemical data
- Sander et al. (2011): Chemical kinetics evaluation

Rate constants are given in the form:
    k = A * (T/300)^n * exp(-Ea/RT)

where:
    A = pre-exponential factor (cm³ molecule⁻¹ s⁻¹ for bimolecular)
    n = temperature exponent
    Ea = activation energy (J/mol)
    R = 8.314 J mol⁻¹ K⁻¹
    T = temperature (K)

Author: GECKO-A Development Team
"""

import math
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum


# ==============================================================================
# Constants
# ==============================================================================

R_GAS = 8.314  # J/(mol·K)
AVOGADRO = 6.022e23  # molecules/mol


# ==============================================================================
# Enums and Data Classes
# ==============================================================================

class OxidantType(Enum):
    """Types of atmospheric oxidants."""
    OH = "OH"  # Hydroxyl radical
    O3 = "O3"  # Ozone
    NO3 = "NO3"  # Nitrate radical
    CL = "Cl"  # Chlorine atom
    O3P = "O(³P)"  # Ground state oxygen atom
    O1D = "O(¹D)"  # Excited state oxygen atom


class ReactionType(Enum):
    """Types of atmospheric reactions."""
    ABSTRACTION = "H-abstraction"
    ADDITION = "addition"
    DECOMPOSITION = "decomposition"
    ISOMERIZATION = "isomerization"
    PHOTOLYSIS = "photolysis"
    ASSOCIATION = "association"


@dataclass
class ArrheniusParams:
    """Arrhenius parameters for rate constant calculation."""
    A: float  # Pre-exponential factor (cm³ molecule⁻¹ s⁻¹)
    n: float = 0.0  # Temperature exponent
    Ea: float = 0.0  # Activation energy (J/mol)
    uncertainty_factor: float = 1.3  # Uncertainty factor (1σ)
    reference: str = ""
    temperature_range: Tuple[float, float] = (200.0, 400.0)  # Valid T range (K)

    def calculate_k(self, T: float = 298.0) -> float:
        """Calculate rate constant at temperature T."""
        return self.A * (T / 300.0) ** self.n * math.exp(-self.Ea / (R_GAS * T))

    def calculate_k_uncertainty(self, T: float = 298.0) -> Tuple[float, float]:
        """Calculate rate constant with uncertainty bounds."""
        k = self.calculate_k(T)
        return (k / self.uncertainty_factor, k * self.uncertainty_factor)


@dataclass
class TroePressureParams:
    """Troe formulation for pressure-dependent reactions."""
    k0_A: float  # Low-pressure limit A
    k0_n: float  # Low-pressure limit T-exponent
    k0_Ea: float = 0.0  # Low-pressure limit Ea
    kinf_A: float = 0.0  # High-pressure limit A
    kinf_n: float = 0.0  # High-pressure limit T-exponent
    kinf_Ea: float = 0.0  # High-pressure limit Ea
    Fc: float = 0.6  # Broadening factor
    reference: str = ""

    def calculate_k(self, T: float = 298.0, M: float = 2.5e19) -> float:
        """
        Calculate rate constant using Troe formulation.

        Args:
            T: Temperature (K)
            M: Total number density (molecules/cm³)
               Default is ~1 atm at 298 K
        """
        k0 = self.k0_A * (T / 300.0) ** self.k0_n * math.exp(-self.k0_Ea / (R_GAS * T))
        kinf = self.kinf_A * (T / 300.0) ** self.kinf_n * math.exp(-self.kinf_Ea / (R_GAS * T))

        if kinf == 0:
            return k0 * M

        Pr = k0 * M / kinf
        log_Fc = math.log10(self.Fc) if self.Fc > 0 else -0.22

        # Broadening factor
        N = 0.75 - 1.27 * log_Fc
        d = 0.14
        log_F = log_Fc / (1 + ((math.log10(Pr) + 0.12) / N) ** 2)
        F = 10 ** log_F

        k = kinf * (Pr / (1 + Pr)) * F
        return k


@dataclass
class BranchingChannel:
    """A single product channel with branching ratio."""
    products: List[str]
    branching_ratio: float  # Fraction (0-1)
    temperature_dependence: Optional[ArrheniusParams] = None
    notes: str = ""


@dataclass
class ReactionKinetics:
    """Complete kinetics data for a reaction."""
    reactant: str
    oxidant: OxidantType
    reaction_type: ReactionType
    rate_params: Union[ArrheniusParams, TroePressureParams]
    channels: List[BranchingChannel] = field(default_factory=list)
    soa_yield: float = 0.0  # Approximate SOA yield
    reference: str = ""
    notes: str = ""


# ==============================================================================
# Rate Constant Database
# ==============================================================================

class ReactionDatabase:
    """
    Database of atmospheric reaction kinetics.

    Contains validated rate constants and branching ratios for:
    - VOC + OH reactions
    - VOC + O3 reactions
    - VOC + NO3 reactions
    - RO2 + NO reactions
    - RO2 + HO2 reactions
    """

    def __init__(self):
        self._reactions: Dict[str, Dict[OxidantType, ReactionKinetics]] = {}
        self._load_database()

    def _load_database(self):
        """Load the reaction database."""
        # =====================================================================
        # ALKANE + OH REACTIONS
        # =====================================================================

        # Methane + OH (JPL 19-5)
        self._add_reaction(ReactionKinetics(
            reactant="methane",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=2.45e-12,
                n=0,
                Ea=7548.0,  # 1775 K * R
                uncertainty_factor=1.1,
                reference="JPL 19-5"
            ),
            channels=[
                BranchingChannel(products=["CH3O2"], branching_ratio=1.0)
            ],
            soa_yield=0.0
        ))

        # Ethane + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="ethane",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=1.52e-17,
                n=2.0,
                Ea=1480.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["C2H5O2"], branching_ratio=1.0)
            ],
            soa_yield=0.0
        ))

        # Propane + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="propane",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=1.65e-17,
                n=2.0,
                Ea=460.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["IC3H7O2"], branching_ratio=0.74,
                                notes="Secondary H abstraction"),
                BranchingChannel(products=["NC3H7O2"], branching_ratio=0.26,
                                notes="Primary H abstraction")
            ],
            soa_yield=0.0
        ))

        # n-Butane + OH (MCM)
        self._add_reaction(ReactionKinetics(
            reactant="n-butane",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=1.36e-17,
                n=2.0,
                Ea=-167.0,
                uncertainty_factor=1.2,
                reference="MCM v3.3.1"
            ),
            channels=[
                BranchingChannel(products=["SC4H9O2"], branching_ratio=0.87),
                BranchingChannel(products=["NC4H9O2"], branching_ratio=0.13)
            ],
            soa_yield=0.01
        ))

        # n-Decane + OH (Atkinson 2003)
        self._add_reaction(ReactionKinetics(
            reactant="n-decane",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=1.10e-11,
                n=0,
                Ea=0,
                uncertainty_factor=1.2,
                reference="Atkinson 2003"
            ),
            soa_yield=0.15
        ))

        # =====================================================================
        # AROMATIC + OH REACTIONS
        # =====================================================================

        # Benzene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="benzene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.33e-12,
                n=0,
                Ea=786.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["phenol", "HO2"], branching_ratio=0.53),
                BranchingChannel(products=["ring-opened products"], branching_ratio=0.47)
            ],
            soa_yield=0.30
        ))

        # Toluene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="toluene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.81e-12,
                n=0,
                Ea=-338.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["cresol", "HO2"], branching_ratio=0.18,
                                notes="Ring addition"),
                BranchingChannel(products=["benzaldehyde", "HO2"], branching_ratio=0.07,
                                notes="Side-chain abstraction"),
                BranchingChannel(products=["ring-opened products"], branching_ratio=0.75)
            ],
            soa_yield=0.35
        ))

        # m-Xylene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="m-xylene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.31e-11,
                n=0,
                Ea=0,
                uncertainty_factor=1.2,
                reference="IUPAC"
            ),
            soa_yield=0.10
        ))

        # =====================================================================
        # MONOTERPENE + OH REACTIONS
        # =====================================================================

        # α-Pinene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="alpha-pinene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.21e-11,
                n=0,
                Ea=-444.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["pinonaldehyde", "OH"], branching_ratio=0.30),
                BranchingChannel(products=["pinene-derived RO2"], branching_ratio=0.70)
            ],
            soa_yield=0.35
        ))

        # β-Pinene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="beta-pinene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.55e-11,
                n=0,
                Ea=-467.0,
                uncertainty_factor=1.2,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["nopinone", "formaldehyde"], branching_ratio=0.35),
                BranchingChannel(products=["beta-pinene RO2"], branching_ratio=0.65)
            ],
            soa_yield=0.30
        ))

        # Limonene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="limonene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.70e-10,
                n=0,
                Ea=0,
                uncertainty_factor=1.3,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["limonaldehyde", "OH"], branching_ratio=0.25),
                BranchingChannel(products=["limonene-derived RO2"], branching_ratio=0.75)
            ],
            soa_yield=0.40
        ))

        # Myrcene + OH (Atkinson)
        self._add_reaction(ReactionKinetics(
            reactant="myrcene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.15e-10,
                n=0,
                Ea=0,
                uncertainty_factor=1.3,
                reference="Atkinson 1997"
            ),
            soa_yield=0.35
        ))

        # =====================================================================
        # MONOTERPENE + O3 REACTIONS
        # =====================================================================

        # α-Pinene + O3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="alpha-pinene",
            oxidant=OxidantType.O3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=6.30e-16,
                n=0,
                Ea=1966.0,
                uncertainty_factor=1.25,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["pinonaldehyde", "OH"], branching_ratio=0.80,
                                notes="OH yield ~0.8"),
                BranchingChannel(products=["stabilized Criegee"], branching_ratio=0.20)
            ],
            soa_yield=0.25
        ))

        # β-Pinene + O3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="beta-pinene",
            oxidant=OxidantType.O3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.50e-17,
                n=0,
                Ea=0,
                uncertainty_factor=1.5,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["nopinone", "formaldehyde", "OH"],
                                branching_ratio=0.35)
            ],
            soa_yield=0.20
        ))

        # Limonene + O3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="limonene",
            oxidant=OxidantType.O3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.80e-15,
                n=0,
                Ea=1280.0,
                uncertainty_factor=1.3,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["limonaldehyde", "OH"], branching_ratio=0.67)
            ],
            soa_yield=0.35
        ))

        # =====================================================================
        # ISOPRENE REACTIONS
        # =====================================================================

        # Isoprene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="isoprene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.70e-11,
                n=0,
                Ea=-390.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["MVK", "formaldehyde", "HO2"],
                                branching_ratio=0.32,
                                notes="4,3-addition"),
                BranchingChannel(products=["MACR", "formaldehyde", "HO2"],
                                branching_ratio=0.23,
                                notes="1,2-addition"),
                BranchingChannel(products=["isoprene-derived RO2"],
                                branching_ratio=0.45)
            ],
            soa_yield=0.03,
            notes="Low SOA yield under high-NOx, higher under low-NOx"
        ))

        # Isoprene + O3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="isoprene",
            oxidant=OxidantType.O3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.03e-14,
                n=0,
                Ea=1995.0,
                uncertainty_factor=1.25,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["MVK", "formaldehyde"], branching_ratio=0.40),
                BranchingChannel(products=["MACR", "formaldehyde"], branching_ratio=0.20),
                BranchingChannel(products=["stabilized Criegee"], branching_ratio=0.40)
            ],
            soa_yield=0.01
        ))

        # Isoprene + NO3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="isoprene",
            oxidant=OxidantType.NO3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=3.03e-12,
                n=0,
                Ea=446.0,
                uncertainty_factor=1.5,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["isoprene nitrates"], branching_ratio=0.80)
            ],
            soa_yield=0.15,
            notes="Nighttime chemistry, organic nitrates are major products"
        ))

        # =====================================================================
        # SESQUITERPENE REACTIONS
        # =====================================================================

        # β-Caryophyllene + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="beta-caryophyllene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.00e-10,
                n=0,
                Ea=0,
                uncertainty_factor=1.5,
                reference="IUPAC"
            ),
            soa_yield=0.60
        ))

        # β-Caryophyllene + O3 (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="beta-caryophyllene",
            oxidant=OxidantType.O3,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=1.16e-14,
                n=0,
                Ea=0,
                uncertainty_factor=1.3,
                reference="IUPAC"
            ),
            soa_yield=0.50
        ))

        # α-Humulene + OH
        self._add_reaction(ReactionKinetics(
            reactant="alpha-humulene",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.93e-10,
                n=0,
                Ea=0,
                uncertainty_factor=1.5,
                reference="Atkinson 1997"
            ),
            soa_yield=0.55
        ))

        # =====================================================================
        # OXYGENATED VOC REACTIONS
        # =====================================================================

        # Formaldehyde + OH (JPL)
        self._add_reaction(ReactionKinetics(
            reactant="formaldehyde",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=5.50e-12,
                n=0,
                Ea=-498.0,
                uncertainty_factor=1.15,
                reference="JPL 19-5"
            ),
            channels=[
                BranchingChannel(products=["HCO", "H2O"], branching_ratio=1.0)
            ],
            soa_yield=0.0
        ))

        # Acetaldehyde + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="acetaldehyde",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=4.40e-12,
                n=0,
                Ea=-365.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["CH3CO", "H2O"], branching_ratio=0.95),
                BranchingChannel(products=["CH2CHO", "H2O"], branching_ratio=0.05)
            ],
            soa_yield=0.0
        ))

        # Methanol + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="methanol",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=2.90e-12,
                n=0,
                Ea=720.0,
                uncertainty_factor=1.15,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["CH2OH", "H2O"], branching_ratio=0.85),
                BranchingChannel(products=["CH3O", "H2O"], branching_ratio=0.15)
            ],
            soa_yield=0.0
        ))

        # Acetone + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="acetone",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ABSTRACTION,
            rate_params=ArrheniusParams(
                A=4.56e-14,
                n=3.65,
                Ea=-2070.0,
                uncertainty_factor=1.2,
                reference="IUPAC"
            ),
            channels=[
                BranchingChannel(products=["CH3COCH2O2"], branching_ratio=1.0)
            ],
            soa_yield=0.0
        ))

        # MVK + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="methyl vinyl ketone",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=2.60e-12,
                n=0,
                Ea=-610.0,
                uncertainty_factor=1.2,
                reference="IUPAC"
            ),
            soa_yield=0.05
        ))

        # MACR + OH (IUPAC)
        self._add_reaction(ReactionKinetics(
            reactant="methacrolein",
            oxidant=OxidantType.OH,
            reaction_type=ReactionType.ADDITION,
            rate_params=ArrheniusParams(
                A=8.00e-12,
                n=0,
                Ea=-380.0,
                uncertainty_factor=1.2,
                reference="IUPAC"
            ),
            soa_yield=0.05
        ))

    def _add_reaction(self, kinetics: ReactionKinetics):
        """Add a reaction to the database."""
        reactant = kinetics.reactant.lower().replace('-', '').replace(' ', '')

        if reactant not in self._reactions:
            self._reactions[reactant] = {}

        self._reactions[reactant][kinetics.oxidant] = kinetics

    def get_reaction(self, reactant: str, oxidant: OxidantType) -> Optional[ReactionKinetics]:
        """Get kinetics for a specific reaction."""
        reactant_key = reactant.lower().replace('-', '').replace(' ', '')

        if reactant_key in self._reactions:
            return self._reactions[reactant_key].get(oxidant)
        return None

    def get_all_reactions(self, reactant: str) -> Dict[OxidantType, ReactionKinetics]:
        """Get all reactions for a given reactant."""
        reactant_key = reactant.lower().replace('-', '').replace(' ', '')
        return self._reactions.get(reactant_key, {})

    def list_reactants(self) -> List[str]:
        """List all reactants in the database."""
        return list(self._reactions.keys())


# ==============================================================================
# Helper Functions
# ==============================================================================

def get_rate_constant(reactant: str, oxidant: str, T: float = 298.0,
                      db: Optional[ReactionDatabase] = None) -> Optional[float]:
    """
    Get rate constant for a reaction at specified temperature.

    Args:
        reactant: Reactant name (e.g., "alpha-pinene", "toluene")
        oxidant: Oxidant name ("OH", "O3", "NO3")
        T: Temperature in Kelvin
        db: Optional ReactionDatabase instance (uses module-level singleton if not provided)

    Returns:
        Rate constant in cm³ molecule⁻¹ s⁻¹, or None if not found
    """
    # Use provided database or defer to module-level singleton (initialized at module load)
    if db is None:
        # Import at function scope to avoid circular import during module initialization
        db = globals().get('reaction_database')
        if db is None:
            db = ReactionDatabase()

    # Map oxidant string to enum
    oxidant_map = {
        'OH': OxidantType.OH,
        'O3': OxidantType.O3,
        'NO3': OxidantType.NO3,
        'Cl': OxidantType.CL,
    }

    if oxidant.upper() not in oxidant_map:
        return None

    kinetics = db.get_reaction(reactant, oxidant_map[oxidant.upper()])
    if kinetics is None:
        return None

    return kinetics.rate_params.calculate_k(T)


def get_branching_ratios(reactant: str, oxidant: str,
                         db: Optional[ReactionDatabase] = None) -> Optional[Dict[str, float]]:
    """
    Get product branching ratios for a reaction.

    Args:
        reactant: Reactant name
        oxidant: Oxidant name
        db: Optional ReactionDatabase instance (uses module-level singleton if not provided)

    Returns:
        Dictionary of products to branching ratios, or None if not found
    """
    if db is None:
        db = globals().get('reaction_database')
        if db is None:
            db = ReactionDatabase()

    oxidant_map = {
        'OH': OxidantType.OH,
        'O3': OxidantType.O3,
        'NO3': OxidantType.NO3,
    }

    if oxidant.upper() not in oxidant_map:
        return None

    kinetics = db.get_reaction(reactant, oxidant_map[oxidant.upper()])
    if kinetics is None or not kinetics.channels:
        return None

    result = {}
    for channel in kinetics.channels:
        products_str = " + ".join(channel.products)
        result[products_str] = channel.branching_ratio

    return result


def estimate_soa_yield(reactant: str, oxidant: str = "OH",
                       db: Optional[ReactionDatabase] = None) -> float:
    """
    Get estimated SOA yield for a reaction.

    Args:
        reactant: Reactant name
        oxidant: Oxidant name (default "OH")
        db: Optional ReactionDatabase instance (uses module-level singleton if not provided)

    Returns:
        Estimated SOA mass yield (0-1)
    """
    if db is None:
        db = globals().get('reaction_database')
        if db is None:
            db = ReactionDatabase()

    oxidant_map = {
        'OH': OxidantType.OH,
        'O3': OxidantType.O3,
        'NO3': OxidantType.NO3,
    }

    if oxidant.upper() not in oxidant_map:
        return 0.0

    kinetics = db.get_reaction(reactant, oxidant_map[oxidant.upper()])
    if kinetics is None:
        # Return default estimates based on compound class
        reactant_lower = reactant.lower()
        if 'pinene' in reactant_lower or 'limonene' in reactant_lower:
            return 0.30
        elif 'caryophyllene' in reactant_lower or 'humulene' in reactant_lower:
            return 0.50
        elif 'toluene' in reactant_lower or 'xylene' in reactant_lower:
            return 0.10
        elif 'benzene' in reactant_lower:
            return 0.30
        elif 'isoprene' in reactant_lower:
            return 0.03
        return 0.05  # Default

    return kinetics.soa_yield


def calculate_atmospheric_lifetime(reactant: str, T: float = 298.0,
                                   oh_conc: float = 1e6,
                                   o3_conc: float = 7.5e11,
                                   no3_conc: float = 2e8,
                                   db: Optional[ReactionDatabase] = None) -> Dict[str, float]:
    """
    Calculate atmospheric lifetimes for a VOC.

    Args:
        reactant: Reactant name
        T: Temperature (K)
        oh_conc: OH concentration (molecules/cm³), default = 1×10⁶
        o3_conc: O3 concentration (molecules/cm³), default = 30 ppb
        no3_conc: NO3 concentration (molecules/cm³), default = 2×10⁸
        db: Optional ReactionDatabase instance (uses module-level singleton if not provided)

    Returns:
        Dictionary of lifetimes in hours for each oxidant
    """
    if db is None:
        db = globals().get('reaction_database')
        if db is None:
            db = ReactionDatabase()
    lifetimes = {}

    oxidants = {
        'OH': (OxidantType.OH, oh_conc),
        'O3': (OxidantType.O3, o3_conc),
        'NO3': (OxidantType.NO3, no3_conc),
    }

    for name, (oxidant, conc) in oxidants.items():
        kinetics = db.get_reaction(reactant, oxidant)
        if kinetics:
            k = kinetics.rate_params.calculate_k(T)
            tau = 1.0 / (k * conc)  # seconds
            lifetimes[name] = tau / 3600.0  # convert to hours
        else:
            lifetimes[name] = float('inf')

    # Calculate total lifetime (1/τ_total = Σ 1/τ_i)
    inv_total = sum(1.0/t for t in lifetimes.values() if t != float('inf'))
    lifetimes['total'] = 1.0/inv_total if inv_total > 0 else float('inf')

    return lifetimes


# ==============================================================================
# Module Test
# ==============================================================================

if __name__ == "__main__":
    print("Reaction Database Summary")
    print("=" * 60)

    db = ReactionDatabase()
    reactants = db.list_reactants()

    print(f"Total reactants: {len(reactants)}")
    print("\nSample rate constants at 298 K:")

    for compound in ['alpha-pinene', 'limonene', 'toluene', 'isoprene']:
        compound_key = compound.lower().replace('-', '')
        if compound_key in reactants:
            for oxidant in [OxidantType.OH, OxidantType.O3, OxidantType.NO3]:
                k = get_rate_constant(compound, oxidant.value)
                if k:
                    print(f"  {compound} + {oxidant.value}: k = {k:.2e} cm³/molecule/s")

    print("\nAtmospheric lifetimes for α-pinene:")
    lifetimes = calculate_atmospheric_lifetime('alpha-pinene')
    for oxidant, tau in lifetimes.items():
        if tau != float('inf'):
            print(f"  τ({oxidant}): {tau:.1f} hours")


# ==============================================================================
# Global Database Instance
# ==============================================================================

# Singleton instance for module-level access
reaction_database = ReactionDatabase()
