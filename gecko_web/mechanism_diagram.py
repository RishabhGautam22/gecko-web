"""
GECKO-A Mechanism Diagram Generator

This module generates comprehensive reaction mechanism outputs including:
1. Static diagrams (PNG/SVG) showing hierarchical oxidation pathways
2. Cytoscape.js-compatible JSON for interactive visualization
3. KPP format for kinetics preprocessor
4. MCM format for Master Chemical Mechanism compatibility
5. FACSIMILE format for atmospheric chemistry modeling

The SMILES conversion uses the unified authoritative function from reaction_tree.py,
which implements proper RDKit validation and an expanded KNOWN_SPECIES database.

Author: Deeksha Sharma
"""

import os
import re
import json
import logging
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from collections import defaultdict

# Import the unified SMILES conversion function from reaction_tree
from .reaction_tree import gecko_to_smiles as unified_gecko_to_smiles, KNOWN_SPECIES as UNIFIED_KNOWN_SPECIES
from .reaction_tree import identify_functional_groups as unified_identify_functional_groups
from .reaction_tree import get_compound_name, GECKO_CODE_TO_NAME

logger = logging.getLogger(__name__)

# Optional dependencies
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available - static diagram generation disabled")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.warning("RDKit not available - will use built-in SMILES generation")


# ==============================================================================
# Data Classes
# ==============================================================================

@dataclass
class Atom:
    """Represents an atom in a molecular graph."""
    index: int
    element: str  # C, O, N, H, S, etc.
    hybridization: str = "sp3"  # sp3, sp2, sp, aromatic
    implicit_h: int = 0
    charge: int = 0
    is_radical: bool = False
    ring_markers: List[int] = field(default_factory=list)


@dataclass
class Bond:
    """Represents a bond between two atoms."""
    atom1_idx: int
    atom2_idx: int
    order: int = 1  # 1=single, 2=double, 3=triple, 4=aromatic


class MolecularGraph:
    """
    Graph representation of a molecule.
    Properly handles branching, ring closures, and bond orders.
    """

    def __init__(self):
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self._adjacency: Dict[int, List[Tuple[int, int]]] = defaultdict(list)  # idx -> [(neighbor_idx, bond_order)]
        self._ring_markers: Dict[int, int] = {}  # marker_num -> atom_idx

    def add_atom(self, element: str, hybridization: str = "sp3",
                 implicit_h: int = 0, is_radical: bool = False) -> int:
        """Add an atom and return its index."""
        idx = len(self.atoms)
        atom = Atom(
            index=idx,
            element=element,
            hybridization=hybridization,
            implicit_h=implicit_h,
            is_radical=is_radical
        )
        self.atoms.append(atom)
        return idx

    def add_bond(self, atom1_idx: int, atom2_idx: int, order: int = 1):
        """Add a bond between two atoms."""
        if atom1_idx < 0 or atom1_idx >= len(self.atoms):
            return
        if atom2_idx < 0 or atom2_idx >= len(self.atoms):
            return
        if atom1_idx == atom2_idx:
            return

        # Check if bond already exists
        for neighbor_idx, _ in self._adjacency[atom1_idx]:
            if neighbor_idx == atom2_idx:
                return

        bond = Bond(atom1_idx, atom2_idx, order)
        self.bonds.append(bond)
        self._adjacency[atom1_idx].append((atom2_idx, order))
        self._adjacency[atom2_idx].append((atom1_idx, order))

    def set_ring_marker(self, marker_num: int, atom_idx: int):
        """Set or resolve a ring marker."""
        if marker_num in self._ring_markers:
            # Close the ring
            other_idx = self._ring_markers[marker_num]
            self.add_bond(atom_idx, other_idx, 1)
            del self._ring_markers[marker_num]
        else:
            self._ring_markers[marker_num] = atom_idx

    def to_smiles(self) -> str:
        """
        Convert molecular graph to SMILES using DFS traversal.
        """
        if not self.atoms:
            return ""

        # Build adjacency matrix with bond orders
        visited = set()
        smiles_parts = []

        def dfs(atom_idx: int, parent_idx: int = -1) -> str:
            if atom_idx in visited:
                return ""

            visited.add(atom_idx)
            atom = self.atoms[atom_idx]

            # Build atom string
            atom_str = self._atom_to_smiles(atom)

            # Get neighbors
            neighbors = [(n_idx, order) for n_idx, order in self._adjacency[atom_idx]
                        if n_idx != parent_idx]

            # Separate visited (ring closures) from unvisited (branches)
            ring_closures = [(n_idx, order) for n_idx, order in neighbors if n_idx in visited]
            branches = [(n_idx, order) for n_idx, order in neighbors if n_idx not in visited]

            # Add ring closure markers
            for rc_idx, rc_order in ring_closures:
                bond_char = "=" if rc_order == 2 else "#" if rc_order == 3 else ""
                atom_str += bond_char + str(rc_idx % 10)

            # Process branches
            if len(branches) == 0:
                return atom_str
            elif len(branches) == 1:
                n_idx, order = branches[0]
                bond_char = "=" if order == 2 else "#" if order == 3 else ""
                return atom_str + bond_char + dfs(n_idx, atom_idx)
            else:
                result = atom_str
                # First branch is main chain
                n_idx, order = branches[0]
                bond_char = "=" if order == 2 else "#" if order == 3 else ""
                main_branch = bond_char + dfs(n_idx, atom_idx)

                # Other branches in parentheses
                for n_idx, order in branches[1:]:
                    bond_char = "=" if order == 2 else "#" if order == 3 else ""
                    result += "(" + bond_char + dfs(n_idx, atom_idx) + ")"

                return result + main_branch

        smiles = dfs(0)
        return smiles

    def _atom_to_smiles(self, atom: Atom) -> str:
        """Convert single atom to SMILES notation."""
        element = atom.element

        # Organic subset: C, N, O, S, P, F, Cl, Br, I
        organic_subset = {'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'}

        if atom.hybridization == "aromatic":
            return element.lower()

        # Need brackets for charged, radical, or non-organic atoms
        needs_bracket = (
            atom.charge != 0 or
            atom.is_radical or
            element not in organic_subset
        )

        if needs_bracket:
            result = "["
            result += element
            if atom.is_radical:
                result += "."
            if atom.charge > 0:
                result += "+" if atom.charge == 1 else f"+{atom.charge}"
            elif atom.charge < 0:
                result += "-" if atom.charge == -1 else f"{atom.charge}"
            result += "]"
            return result

        return element


# ==============================================================================
# Graph-Based GECKO Formula Parser
# ==============================================================================

class GeckoGraphParser:
    """
    Parses GECKO-A structural formulas into MolecularGraph objects.

    This class now uses the unified SMILES conversion from reaction_tree.py
    which has an expanded KNOWN_SPECIES database and proper RDKit validation.
    """

    # Use the unified KNOWN_SPECIES database from reaction_tree
    KNOWN_SPECIES = UNIFIED_KNOWN_SPECIES

    @classmethod
    def to_smiles(cls, gecko_formula: str, species_code: str = "") -> str:
        """
        Convert GECKO-A formula to valid SMILES.

        This method now delegates to the unified gecko_to_smiles function
        from reaction_tree.py for consistent SMILES generation across the application.

        Args:
            gecko_formula: The GECKO formula string
            species_code: Optional GECKO species code

        Returns:
            Valid SMILES string
        """
        # Delegate to the unified function
        return unified_gecko_to_smiles(gecko_formula, species_code)

    @classmethod
    def _parse_to_graph(cls, formula: str) -> MolecularGraph:
        """
        Parse GECKO formula into a molecular graph using stack-based branch tracking.
        """
        graph = MolecularGraph()
        formula = formula.replace('·', '.').replace('•', '.')

        # State tracking
        current_atom_idx = -1  # Index of current atom to attach to
        branch_stack = []  # Stack of attachment points for branching
        pending_double_bond = False  # Flag for Cd-Cd double bonds
        i = 0
        n = len(formula)

        while i < n:
            char = formula[i]

            # Skip single dashes (bond separators)
            if char == '-' and (i + 1 >= n or formula[i + 1] != '-'):
                i += 1
                continue

            # Double bond marker
            if char == '=':
                pending_double_bond = True
                i += 1
                continue

            # Opening parenthesis - start branch or functional group
            if char == '(':
                # Find matching close paren
                depth = 1
                j = i + 1
                while j < n and depth > 0:
                    if formula[j] == '(':
                        depth += 1
                    elif formula[j] == ')':
                        depth -= 1
                    j += 1

                group_content = formula[i + 1:j - 1]
                group_upper = group_content.upper().strip()

                # Check if it's a known functional group
                func_atoms = cls._parse_functional_group(group_upper, graph, current_atom_idx)
                if func_atoms is not None:
                    i = j
                    continue

                # It's a branch - push current position and parse recursively
                if current_atom_idx >= 0:
                    branch_stack.append(current_atom_idx)

                # Parse branch content
                branch_graph = cls._parse_to_graph(group_content)
                if branch_graph and branch_graph.atoms:
                    # Merge branch into main graph
                    offset = len(graph.atoms)
                    for atom in branch_graph.atoms:
                        new_idx = graph.add_atom(
                            atom.element, atom.hybridization,
                            atom.implicit_h, atom.is_radical
                        )

                    # Copy bonds with offset
                    for bond in branch_graph.bonds:
                        graph.add_bond(bond.atom1_idx + offset, bond.atom2_idx + offset, bond.order)

                    # Connect first atom of branch to current atom
                    if current_atom_idx >= 0:
                        bond_order = 2 if pending_double_bond else 1
                        graph.add_bond(current_atom_idx, offset, bond_order)
                        pending_double_bond = False

                # Pop branch stack (stay at same current_atom for next branch or main chain)
                if branch_stack:
                    current_atom_idx = branch_stack.pop()

                i = j
                continue

            # Closing parenthesis (shouldn't happen in balanced input, but handle gracefully)
            if char == ')':
                i += 1
                continue

            # Carbon groups
            if char == 'C':
                # Check for Cd (sp2 carbon)
                if i + 1 < n and formula[i + 1] == 'd':
                    hybridization = "sp2"
                    i += 2

                    # Check for CdH, CdH2, CdH3
                    implicit_h = 0
                    if i < n and formula[i] == 'H':
                        i += 1
                        if i < n and formula[i].isdigit():
                            implicit_h = int(formula[i])
                            i += 1
                        else:
                            implicit_h = 1

                    new_idx = graph.add_atom('C', hybridization, implicit_h)

                    if current_atom_idx >= 0:
                        # If both are sp2 and no explicit bond, make double bond
                        bond_order = 1
                        if pending_double_bond:
                            bond_order = 2
                            pending_double_bond = False
                        elif graph.atoms[current_atom_idx].hybridization == "sp2":
                            bond_order = 2

                        graph.add_bond(current_atom_idx, new_idx, bond_order)

                    current_atom_idx = new_idx
                    continue

                # Check for ring marker (C1, C2, etc.)
                if i + 1 < n and formula[i + 1].isdigit():
                    ring_num = int(formula[i + 1])
                    new_idx = graph.add_atom('C', 'sp3', 0)

                    if current_atom_idx >= 0:
                        bond_order = 2 if pending_double_bond else 1
                        graph.add_bond(current_atom_idx, new_idx, bond_order)
                        pending_double_bond = False

                    graph.set_ring_marker(ring_num, new_idx)
                    current_atom_idx = new_idx
                    i += 2
                    continue

                # Check for CHn
                if i + 1 < n and formula[i + 1] == 'H':
                    i += 2
                    implicit_h = 1
                    if i < n and formula[i].isdigit():
                        implicit_h = int(formula[i])
                        i += 1

                    new_idx = graph.add_atom('C', 'sp3', implicit_h)

                    if current_atom_idx >= 0:
                        bond_order = 2 if pending_double_bond else 1
                        graph.add_bond(current_atom_idx, new_idx, bond_order)
                        pending_double_bond = False

                    current_atom_idx = new_idx
                    continue

                # Check for CO (carbonyl in chain)
                if i + 1 < n and formula[i + 1] == 'O':
                    # Differentiate C-O-ring from C=O
                    if i + 2 < n and formula[i + 2].isdigit():
                        # It's a ring C-O marker
                        pass
                    else:
                        # Carbonyl group
                        c_idx = graph.add_atom('C', 'sp2', 0)
                        o_idx = graph.add_atom('O', 'sp2', 0)
                        graph.add_bond(c_idx, o_idx, 2)

                        if current_atom_idx >= 0:
                            bond_order = 2 if pending_double_bond else 1
                            graph.add_bond(current_atom_idx, c_idx, bond_order)
                            pending_double_bond = False

                        current_atom_idx = c_idx
                        i += 2
                        continue

                # Plain C
                new_idx = graph.add_atom('C', 'sp3', 0)
                if current_atom_idx >= 0:
                    bond_order = 2 if pending_double_bond else 1
                    graph.add_bond(current_atom_idx, new_idx, bond_order)
                    pending_double_bond = False

                current_atom_idx = new_idx
                i += 1
                continue

            # Oxygen (with potential ring marker)
            if char == 'O':
                if i + 1 < n and formula[i + 1].isdigit():
                    ring_num = int(formula[i + 1])
                    new_idx = graph.add_atom('O', 'sp3', 0)

                    if current_atom_idx >= 0:
                        graph.add_bond(current_atom_idx, new_idx, 1)

                    graph.set_ring_marker(ring_num, new_idx)
                    current_atom_idx = new_idx
                    i += 2
                    continue

                # Plain O
                new_idx = graph.add_atom('O', 'sp3', 0)
                if current_atom_idx >= 0:
                    bond_order = 2 if pending_double_bond else 1
                    graph.add_bond(current_atom_idx, new_idx, bond_order)
                    pending_double_bond = False

                current_atom_idx = new_idx
                i += 1
                continue

            # Aromatic carbon (lowercase c)
            if char == 'c':
                new_idx = graph.add_atom('C', 'aromatic', 0)

                if current_atom_idx >= 0:
                    graph.add_bond(current_atom_idx, new_idx, 4)  # 4 = aromatic

                # Check for ring marker
                if i + 1 < n and formula[i + 1].isdigit():
                    ring_num = int(formula[i + 1])
                    graph.set_ring_marker(ring_num, new_idx)
                    i += 2
                else:
                    i += 1

                # Skip attached H
                if i < n and formula[i] == 'H':
                    i += 1

                current_atom_idx = new_idx
                continue

            # Nitrogen
            if char == 'N':
                new_idx = graph.add_atom('N', 'sp3', 0)
                if current_atom_idx >= 0:
                    bond_order = 2 if pending_double_bond else 1
                    graph.add_bond(current_atom_idx, new_idx, bond_order)
                    pending_double_bond = False

                current_atom_idx = new_idx
                i += 1
                continue

            # Skip unknown characters
            i += 1

        return graph

    @classmethod
    def _parse_functional_group(cls, group: str, graph: MolecularGraph,
                                attach_to: int) -> Optional[int]:
        """
        Parse a known functional group and attach it to the graph.
        Returns the index of the attachment atom, or None if not a known group.
        """
        # Mapping of functional groups to (atoms, bonds) structure
        # Returns index of last atom added for chaining

        if group == 'OH':
            o_idx = graph.add_atom('O', 'sp3', 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, o_idx, 1)
            return o_idx

        if group == 'OOH':
            o1_idx = graph.add_atom('O', 'sp3', 0)
            o2_idx = graph.add_atom('O', 'sp3', 1)
            graph.add_bond(o1_idx, o2_idx, 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, o1_idx, 1)
            return o2_idx

        if group in ('ONO2', 'NO3'):
            o1_idx = graph.add_atom('O', 'sp3', 0)
            n_idx = graph.add_atom('N', 'sp2', 0)
            n_idx_atom = graph.atoms[n_idx]
            n_idx_atom.charge = 1
            o2_idx = graph.add_atom('O', 'sp2', 0)
            o3_idx = graph.add_atom('O', 'sp2', 0)
            o3_atom = graph.atoms[o3_idx]
            o3_atom.charge = -1
            graph.add_bond(o1_idx, n_idx, 1)
            graph.add_bond(n_idx, o2_idx, 2)
            graph.add_bond(n_idx, o3_idx, 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, o1_idx, 1)
            return o1_idx

        if group in ('OO.', 'OO'):
            o1_idx = graph.add_atom('O', 'sp3', 0)
            o2_idx = graph.add_atom('O', 'sp3', 0, is_radical=(group == 'OO.'))
            graph.add_bond(o1_idx, o2_idx, 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, o1_idx, 1)
            return o2_idx

        if group == 'O.':
            o_idx = graph.add_atom('O', 'sp3', 0, is_radical=True)
            if attach_to >= 0:
                graph.add_bond(attach_to, o_idx, 1)
            return o_idx

        if group == 'CO':
            c_idx = graph.add_atom('C', 'sp2', 0)
            o_idx = graph.add_atom('O', 'sp2', 0)
            graph.add_bond(c_idx, o_idx, 2)
            if attach_to >= 0:
                graph.add_bond(attach_to, c_idx, 1)
            return c_idx

        if group == 'CHO':
            c_idx = graph.add_atom('C', 'sp2', 1)
            o_idx = graph.add_atom('O', 'sp2', 0)
            graph.add_bond(c_idx, o_idx, 2)
            if attach_to >= 0:
                graph.add_bond(attach_to, c_idx, 1)
            return c_idx

        if group == 'COOH':
            c_idx = graph.add_atom('C', 'sp2', 0)
            o1_idx = graph.add_atom('O', 'sp2', 0)
            o2_idx = graph.add_atom('O', 'sp3', 1)
            graph.add_bond(c_idx, o1_idx, 2)
            graph.add_bond(c_idx, o2_idx, 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, c_idx, 1)
            return c_idx

        if group == 'NO2':
            n_idx = graph.add_atom('N', 'sp2', 0)
            n_atom = graph.atoms[n_idx]
            n_atom.charge = 1
            o1_idx = graph.add_atom('O', 'sp2', 0)
            o2_idx = graph.add_atom('O', 'sp2', 0)
            o2_atom = graph.atoms[o2_idx]
            o2_atom.charge = -1
            graph.add_bond(n_idx, o1_idx, 2)
            graph.add_bond(n_idx, o2_idx, 1)
            if attach_to >= 0:
                graph.add_bond(attach_to, n_idx, 1)
            return n_idx

        # Handle simple CH3, CH2, CH as branches
        if group in ('CH3', 'CH2', 'CH'):
            h_count = {'CH3': 3, 'CH2': 2, 'CH': 1}.get(group, 0)
            c_idx = graph.add_atom('C', 'sp3', h_count)
            if attach_to >= 0:
                graph.add_bond(attach_to, c_idx, 1)
            return c_idx

        return None

    @classmethod
    def _generate_fallback_structure(cls, formula: str) -> str:
        """Generate a representative SMILES when parsing fails."""
        # Count approximate carbons
        n_carbon = len(re.findall(r'C[dH123]*', formula))
        if n_carbon == 0:
            match = re.search(r'C(\d+)', formula)
            if match:
                n_carbon = int(match.group(1))

        n_carbon = max(1, min(n_carbon, 15))

        # Detect features
        has_double = 'Cd' in formula or '=C' in formula
        has_ring = bool(re.search(r'[CO]\d.*[CO]\d', formula))
        has_aromatic = 'cH' in formula

        has_oh = '(OH)' in formula
        has_ooh = '(OOH)' in formula
        has_ono2 = '(ONO2)' in formula
        has_peroxy = '(OO.)' in formula or 'OO.' in formula
        has_alkoxy = '(O.)' in formula and not has_peroxy
        has_carbonyl = '(CO)' in formula or '-CO-' in formula
        has_aldehyde = 'CHO' in formula
        has_cooh = '(COOH)' in formula

        # Build base structure
        if has_aromatic and n_carbon >= 6:
            base = 'c1ccccc1'
            if n_carbon > 6:
                base = 'C' * (n_carbon - 6) + base
        elif has_ring and n_carbon >= 4:
            ring_size = min(6, n_carbon)
            base = 'C1' + 'C' * (ring_size - 2) + 'C1'
            if n_carbon > ring_size:
                base += 'C' * (n_carbon - ring_size)
        elif has_double and n_carbon >= 2:
            if n_carbon == 5:
                base = 'CC(=C)C=C'  # Isoprene-like
            elif n_carbon >= 4:
                half = n_carbon // 2
                base = 'C' * half + '=C' + 'C' * (n_carbon - half - 1)
            else:
                base = 'C=C'
        else:
            base = 'C' * n_carbon

        # Add functional groups
        suffix = ''
        if has_oh:
            suffix += 'O'
        if has_ooh:
            suffix += 'OO'
        if has_ono2:
            suffix += 'O[N+](=O)[O-]'
        if has_peroxy:
            suffix += 'O[O]'
        if has_alkoxy:
            suffix += '[O]'
        if has_cooh:
            suffix += 'C(=O)O'

        if has_carbonyl or has_aldehyde:
            if len(base) >= 3:
                mid = len(base) // 2
                base = base[:mid] + 'C(=O)' + base[mid + 1:]
            else:
                base = 'CC(=O)C'[:n_carbon * 2]

        return base + suffix

    @classmethod
    def identify_functional_groups(cls, gecko_formula: str) -> List[str]:
        """Identify functional groups in GECKO formula."""
        groups = []
        checks = [
            ('(OH)', 'hydroxyl'),
            ('(OOH)', 'hydroperoxide'),
            ('(ONO2)', 'nitrate'),
            ('(OO.)', 'peroxy radical'),
            ('OO.', 'peroxy radical'),
            ('(O.)', 'alkoxy radical'),
            ('(CO)', 'carbonyl'),
            ('CHO', 'aldehyde'),
            ('(COOH)', 'carboxylic acid'),
            ('(NO2)', 'nitro'),
            ('PAN', 'peroxyacyl nitrate'),
        ]
        for pattern, name in checks:
            if pattern in gecko_formula and name not in groups:
                groups.append(name)
        return groups

    @classmethod
    def is_radical(cls, gecko_formula: str) -> bool:
        """Check if species is a radical."""
        return any(ind in gecko_formula for ind in ['.)', '.]', 'O.', 'OO.', '(O.)', '(OO.)'])


# ==============================================================================
# Data Classes for Mechanism
# ==============================================================================

@dataclass
class Species:
    """Represents a chemical species in the mechanism."""
    code: str
    formula: str = ""
    smiles: str = ""
    molecular_weight: float = 0.0
    n_carbon: int = 0
    n_hydrogen: int = 0
    n_oxygen: int = 0
    n_nitrogen: int = 0
    phase: str = "gas"
    is_radical: bool = False
    functional_groups: List[str] = field(default_factory=list)


@dataclass
class Reaction:
    """Represents a chemical reaction."""
    id: int
    reactants: List[Tuple[str, float]]
    products: List[Tuple[str, float]]
    rate_constant_A: float = 0.0
    rate_constant_n: float = 0.0
    rate_constant_Ea: float = 0.0
    reaction_type: str = "standard"
    comment: str = ""

    def get_rate_expression(self, format_type: str = "kpp") -> str:
        """Return rate expression in specified format."""
        if self.reaction_type == "photolysis":
            return f"J({self.id})"

        if format_type == "kpp":
            if self.rate_constant_n == 0 and self.rate_constant_Ea == 0:
                return f"{self.rate_constant_A:.3E}"
            elif self.rate_constant_n == 0:
                return f"{self.rate_constant_A:.3E}*EXP({-self.rate_constant_Ea}/TEMP)"
            else:
                return f"{self.rate_constant_A:.3E}*(TEMP/300)^{self.rate_constant_n}*EXP({-self.rate_constant_Ea}/TEMP)"
        else:
            return f"k = {self.rate_constant_A:.3E}"


@dataclass
class ReactionPathway:
    """Represents a pathway from parent to products."""
    parent: str
    products: List[str]
    branching_ratio: float
    oxidant: str = ""
    reaction_type: str = "oxidation"
    reaction_ids: List[int] = field(default_factory=list)


@dataclass
class MechanismTree:
    """Hierarchical tree structure of the reaction mechanism."""
    root: str
    root_formula: str
    pathways: Dict[str, List[ReactionPathway]] = field(default_factory=dict)
    species_info: Dict[str, Species] = field(default_factory=dict)
    all_reactions: List[Reaction] = field(default_factory=list)


# ==============================================================================
# Mechanism Parser
# ==============================================================================

class MechanismParser:
    """Parser for GECKO-A mechanism files."""

    OXIDANTS = {'GHO', 'GO3', 'GNO3', 'GHO2', 'GNO', 'GNO2', 'HO', 'O3', 'NO3', 'HO2', 'NO', 'NO2'}

    FILTER_SPECIES = {
        'NOTHING', 'TBODY', 'OXYGEN', 'M', 'N2', 'O2', 'H2O',
        'GH2O', 'GO2', 'GN2', 'EXTRA', 'HV', 'FALLOFF',
        'PERO1', 'PERO2', 'PERO3', 'PERO4', 'PERO5', 'PERO6',
        'PERO7', 'PERO8', 'PERO9', 'MEPERO', 'AIN', 'AOU'
    }

    INORGANIC = {
        'GHO', 'GHO2', 'GNO', 'GNO2', 'GNO3', 'GO3', 'GCO', 'GCO2',
        'GH2O', 'GH2O2', 'GO2', 'GSO2', 'GSULF', 'GH2', 'GO3P', 'GO1D',
        'GHNO2', 'GHNO3', 'GHNO4', 'GN2O5', 'GCH4', 'GCH3O2', 'GCH3O',
        'GCH3OOH', 'GCH2O', 'GCH3OH', 'GHCOOH', 'GHCOO2H', 'GHCOO2',
        'HO', 'HO2', 'NO', 'NO2', 'NO3', 'O3', 'CO', 'CO2',
        'H2O', 'H2O2', 'O2', 'SO2', 'H2', 'HNO2', 'HNO3', 'HNO4', 'N2O5',
        'CH4', 'CH3O2', 'CH3O', 'CH3OOH', 'CH2O', 'CH3OH', 'HCOOH'
    }

    # Mapping of common VOC names to GECKO codes
    VOC_NAME_TO_CODE = {
        'alpha-pinene': 'APINEN',
        'alphapinene': 'APINEN',
        'a-pinene': 'APINEN',
        'apinene': 'APINEN',
        'beta-pinene': 'BPINEN',
        'betapinene': 'BPINEN',
        'b-pinene': 'BPINEN',
        'bpinene': 'BPINEN',
        'isoprene': 'ISOPRE',
        'limonene': 'LIMONE',
        'myrcene': 'MYRCEN',
        'ocimene': 'OCIMEN',
        'toluene': 'TOLUEN',
        'benzene': 'BENZEN',
        'ethane': 'C2H6',
        'propane': 'C3H8',
        'butane': 'NC4H10',
        'pentane': 'NC5H12',
        'hexane': 'NC6H14',
        'heptane': 'NC7H16',
        'octane': 'NC8H18',
        'decane': 'NC10H22',
        'dodecane': 'NC12H26',
        'ethene': 'C2H4',
        'propene': 'C3H6',
        'o-xylene': 'OXYL',
        'm-xylene': 'MXYL',
        'p-xylene': 'PXYL',
        'ethylbenzene': 'EBENZ',
    }

    def __init__(self, output_dir: str, voc_name: str):
        self.output_dir = output_dir
        self.voc_name = voc_name
        self.species_dict: Dict[str, Species] = {}
        self.reactions: List[Reaction] = []
        self.reaction_index: Dict[str, List[int]] = defaultdict(list)
        self.code_map: Dict[str, str] = {}

    def parse(self) -> MechanismTree:
        """Parse mechanism files and build tree structure."""
        self._parse_dictionary()
        self._parse_reactions()

        root_code = self._find_root_voc()
        root_formula = self.species_dict.get(root_code, Species(code=root_code)).formula

        pathways = self._build_pathways(root_code)

        return MechanismTree(
            root=root_code,
            root_formula=root_formula,
            pathways=pathways,
            species_info=self.species_dict,
            all_reactions=self.reactions
        )

    def _parse_dictionary(self):
        """Parse dictionary.out file."""
        dict_path = os.path.join(self.output_dir, "dictionary.out")
        if not os.path.exists(dict_path):
            logger.warning(f"Dictionary file not found: {dict_path}")
            return

        with open(dict_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('!') or 'number of record' in line:
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

                code = parts[0]
                formula = parts[1]

                mw, nC, nH, nN, nO = 0.0, 0, 0, 0, 0
                try:
                    for idx, p in enumerate(parts[2:], start=2):
                        try:
                            val = float(p)
                            if '.' in p:
                                mw = val
                                if idx + 2 < len(parts):
                                    nC = int(parts[idx + 2]) if idx + 2 < len(parts) else 0
                                    nH = int(parts[idx + 3]) if idx + 3 < len(parts) else 0
                                    nN = int(parts[idx + 4]) if idx + 4 < len(parts) else 0
                                    nO = int(parts[idx + 5]) if idx + 5 < len(parts) else 0
                                break
                        except ValueError:
                            continue
                except Exception:
                    pass

                smiles = GeckoGraphParser.to_smiles(formula, code)

                # --- FIX FOR ISOPENTANE/C05000 AMBIGUITY ---
                # If we detect that this job is running 'isopentane' and the code is 'C05000'
                # (or any generic C5 code), force the isopentane SMILES.
                
                is_isopentane_job = False
                if self.voc_name:
                    voc_clean = self.voc_name.lower().strip().replace('-', '').replace('_', '')
                    if voc_clean in ['isopentane', '2methylbutane', 'ipentane']:
                        is_isopentane_job = True

                # Added checking for H05000 and GH05000 explicitly as parent candidates if formula is generic
                if is_isopentane_job: 
                     if (code in {'C05000', 'G05000', '205000', 'GH05000', 'H05000'}) and (smiles in {'CCCCC', 'CCC(C)C'} or not smiles):
                        # Force override
                        logger.info(f"Overriding SMILES for {code} in '{self.voc_name}' job: {smiles} -> CC(C)CC")
                        smiles = 'CC(C)CC'
                # -----------------------------------------------

                species = Species(
                    code=code,
                    formula=formula,
                    smiles=smiles,

                    molecular_weight=mw,
                    n_carbon=nC,
                    n_hydrogen=nH,
                    n_nitrogen=nN,
                    n_oxygen=nO,
                    is_radical=GeckoGraphParser.is_radical(formula),
                    functional_groups=unified_identify_functional_groups(smiles, formula)
                )

                self.species_dict[code] = species
                # Register aliases for Gas (G) and Aerosol (A) phases
                # This ensures AH05000 uses the same structure as H05000
                self.species_dict[f'G{code}'] = species
                self.species_dict[f'A{code}'] = species
                self.species_dict[f'AH{code}'] = species # Cover double prefix cases if any
                
                self.code_map[f'G{code}'] = code
                self.code_map[f'A{code}'] = code
                self.code_map[code] = code

    def _parse_reactions(self):
        """Parse reaction mechanism files."""
        import glob

        reaction_files = [
            os.path.join(self.output_dir, "reactions.txt"),
            os.path.join(self.output_dir, "reactions.dum"),
        ]
        reaction_files.extend(glob.glob(os.path.join(self.output_dir, "*.mech")))

        lines = []
        for rf in reaction_files:
            if os.path.exists(rf):
                with open(rf, 'r') as f:
                    lines.extend(f.readlines())

        if not lines:
            logger.warning("No reaction files found")
            return

        reaction_id = 0
        for line in lines:
            if line.strip().startswith('!'):
                continue
            if any(kw in line for kw in ['SIZE', 'SPECIES', 'PHASE:', 'REACTIONS']):
                continue

            reaction = self._parse_reaction_line(line, reaction_id)
            if reaction:
                self.reactions.append(reaction)
                for species, _ in reaction.reactants:
                    self.reaction_index[species].append(reaction_id)
                reaction_id += 1

    def _parse_reaction_line(self, line: str, reaction_id: int) -> Optional[Reaction]:
        """Parse a single reaction line."""
        if '!' in line:
            comment_idx = line.index('!')
            comment = line[comment_idx + 1:].strip()
            line = line[:comment_idx]
        else:
            comment = ""

        line = line.strip()
        if not line:
            return None

        sep = None
        if '=>' in line:
            sep = '=>'
        elif '=' in line and not any(kw in line for kw in ['FALLOFF', 'EXTRA', 'HV']):
            sep = '='

        if not sep:
            return None

        try:
            parts = line.split(sep, 1)
            lhs = parts[0].strip()
            rhs = parts[1].strip()

            reactants = self._parse_species_list(lhs)
            products, rate_A, rate_n, rate_Ea = self._parse_products_and_rate(rhs)

            if not reactants:
                return None

            reaction_type = "standard"
            if 'HV' in line:
                reaction_type = "photolysis"
            elif 'FALLOFF' in line:
                reaction_type = "falloff"
            elif 'EXTRA' in line:
                reaction_type = "special"

            return Reaction(
                id=reaction_id,
                reactants=reactants,
                products=products,
                rate_constant_A=rate_A,
                rate_constant_n=rate_n,
                rate_constant_Ea=rate_Ea,
                reaction_type=reaction_type,
                comment=comment
            )
        except Exception as e:
            logger.debug(f"Failed to parse reaction: {line[:50]}... Error: {e}")
            return None

    def _normalize_species_code(self, code: str) -> str:
        """Normalize species code by stripping G prefix for gas-phase species."""
        if not code:
            return code
        # Strip G prefix but preserve certain codes that start with G
        if code.startswith('G') and len(code) > 1:
            remainder = code[1:]
            # Preserve GECKO special codes like GLYOX, GLYOXAL
            if code.upper() in {'GLYOX', 'GLYOXAL'}:
                return code
            # Strip G prefix if remainder starts with uppercase letter or digit
            # (e.g., GAPINEN -> APINEN, G2T0000 -> 2T0000, GHO -> HO)
            if remainder[0].isupper() or remainder[0].isdigit():
                return remainder
        return code

    def _parse_species_list(self, text: str) -> List[Tuple[str, float]]:
        """Parse a list of species with stoichiometry."""
        result = []
        parts = text.split('+')

        for part in parts:
            part = part.strip()
            if not part:
                continue

            match = re.match(r'^([\d.]+)?\s*(\w+)', part)
            if match:
                coef = float(match.group(1)) if match.group(1) else 1.0
                species = match.group(2)

                if species not in self.FILTER_SPECIES:
                    # Normalize by stripping G prefix
                    normalized = self._normalize_species_code(species)
                    if normalized not in self.FILTER_SPECIES:
                        result.append((normalized, coef))

        return result

    def _parse_products_and_rate(self, text: str) -> Tuple[List[Tuple[str, float]], float, float, float]:
        """Parse products and rate constants from RHS."""
        tokens = text.split()

        products = []
        rate_A, rate_n, rate_Ea = 1.0, 0.0, 0.0

        rate_start = len(tokens)
        num_count = 0
        for i in range(len(tokens) - 1, -1, -1):
            try:
                float(tokens[i].replace('E', 'e').replace('+', '').replace('-', ''))
                num_count += 1
                if num_count >= 3:
                    rate_start = i
                    break
            except ValueError:
                if num_count > 0:
                    rate_start = i + 1
                    break

        products_text = ' '.join(tokens[:rate_start])
        for part in products_text.split('+'):
            part = part.strip()
            if not part:
                continue

            match = re.match(r'^([\d.]+)?\s*(\w+)', part)
            if match:
                coef = float(match.group(1)) if match.group(1) else 1.0
                species = match.group(2)
                if species not in self.FILTER_SPECIES:
                    # Normalize by stripping G prefix
                    normalized = self._normalize_species_code(species)
                    if normalized not in self.FILTER_SPECIES:
                        products.append((normalized, coef))

        rate_tokens = tokens[rate_start:]
        try:
            if len(rate_tokens) >= 1:
                rate_A = float(rate_tokens[0])
            if len(rate_tokens) >= 2:
                rate_n = float(rate_tokens[1])
            if len(rate_tokens) >= 3:
                rate_Ea = float(rate_tokens[2])
        except ValueError:
            pass

        return products, rate_A, rate_n, rate_Ea

    def _find_root_voc(self) -> str:
        """Find the root VOC species.

        Priority order:
        1. Read from listprimary.dat (most reliable - GECKO-A explicitly specifies this)
        2. Check VOC_NAME_TO_CODE mapping
        3. Match against species dictionary by formula or code
        4. Fuzzy matching fallbacks
        """
        voc_lower = self.voc_name.lower().strip()
        voc_clean = voc_lower.replace('-', '').replace('_', '').replace(' ', '')

        # 1. PRIORITY: Read from listprimary.dat - this is the authoritative source
        # GECKO-A writes the primary VOC code here during mechanism generation
        listprimary_path = os.path.join(self.output_dir, "listprimary.dat")
        if os.path.exists(listprimary_path):
            try:
                with open(listprimary_path, 'r') as f:
                    content = f.read().strip()
                    # Format is: " CODE FORMULA" (space-separated)
                    parts = content.split()
                    if parts:
                        primary_code = parts[0]
                        logger.info(f"Found primary VOC from listprimary.dat: {primary_code}")
                        # Verify it exists in our species dict (with or without G prefix)
                        if primary_code in self.species_dict:
                            return primary_code
                        if f'G{primary_code}' in self.species_dict:
                            return primary_code
                        # Even if not in species_dict, trust listprimary.dat
                        return primary_code
            except Exception as e:
                logger.warning(f"Failed to read listprimary.dat: {e}")

        # 2. Check the VOC name to code mapping
        if voc_lower in self.VOC_NAME_TO_CODE:
            code = self.VOC_NAME_TO_CODE[voc_lower]
            if code in self.species_dict:
                return code
            # Also check with G prefix stripped in pathways
            if f'G{code}' in self.species_dict:
                return code

        # 3. Check for exact code match (with and without G prefix)
        for code in self.species_dict:
            code_clean = code.upper().lstrip('G')
            if code_clean == voc_clean.upper():
                return code.lstrip('G') if code.startswith('G') else code

        # 4. Match against formula
        for code, species in self.species_dict.items():
            if species.formula.lower() == voc_lower:
                return code

        # 5. Fuzzy match on code names
        for code in self.species_dict:
            c_upper = code.upper().lstrip('G')
            v_upper = self.voc_name.upper().replace('-', '').replace('_', '')
            if len(c_upper) >= 4 and len(v_upper) >= 4:
                if v_upper.startswith(c_upper) or c_upper.startswith(v_upper):
                    return code.lstrip('G') if code.startswith('G') else code

        # 6. Special fallback for common VOCs
        common_codes = ['APINEN', 'BPINEN', 'ISOPRE', 'LIMONE', 'MYRCEN', 'TOLUEN']
        for code in common_codes:
            if code in self.species_dict:
                if voc_clean in code.lower() or code.lower() in voc_clean:
                    return code

        logger.warning(f"Could not find root VOC for '{self.voc_name}', using as-is")
        return self.voc_name.upper()

    def _is_inorganic(self, code: str) -> bool:
        """Check if a species code is inorganic (should be filtered from pathways)."""
        # Check both with and without G prefix
        return code in self.INORGANIC or f'G{code}' in self.INORGANIC or code.lstrip('G') in self.INORGANIC

    def _build_pathways(self, root_code: str) -> Dict[str, List[ReactionPathway]]:
        """Build reaction pathways from root."""
        pathways: Dict[str, List[ReactionPathway]] = defaultdict(list)

        for reaction in self.reactions:
            for reactant_raw, _ in reaction.reactants:
                # Normalize reactant code
                reactant = self._normalize_species_code(reactant_raw)

                # Skip inorganic and filter species as reactants
                if self._is_inorganic(reactant) or reactant in self.FILTER_SPECIES:
                    continue

                # Filter products
                products = []
                for p, coef in reaction.products:
                    p_normalized = self._normalize_species_code(p)
                    if not self._is_inorganic(p_normalized) and p_normalized not in self.FILTER_SPECIES:
                        products.append(p_normalized)

                if products:
                    oxidant = ""
                    for r, _ in reaction.reactants:
                        if r in self.OXIDANTS or self._normalize_species_code(r) in {'HO', 'O3', 'NO3', 'HO2'}:
                            oxidant = self._normalize_species_code(r)
                            break

                    pathway = ReactionPathway(
                        parent=reactant,
                        products=products,
                        branching_ratio=1.0,
                        oxidant=oxidant,
                        reaction_ids=[reaction.id]
                    )
                    pathways[reactant].append(pathway)

        return dict(pathways)


# ==============================================================================
# Output Generators
# ==============================================================================

def generate_cytoscape_json(tree: MechanismTree, output_path: str,
                           max_depth: int = 3, max_nodes: int = 50,
                           min_branching: float = 0.05):
    """Generate JSON for interactive visualization matching frontend expectations."""
    nodes = []
    edges = []
    visited = set()
    node_ids = set()  # Track which node IDs have been added

    def add_node(code: str, depth: int):
        if code in visited or len(nodes) >= max_nodes:
            return
        if depth > max_depth:
            return

        visited.add(code)
        node_ids.add(code)  # Track the node ID
        species = tree.species_info.get(code, Species(code=code))

        # Clean SMILES for SmilesDrawer compatibility
        clean_smiles = species.smiles
        if clean_smiles:
            # Replace radical dot notations that SmilesDrawer doesn't understand
            clean_smiles = clean_smiles.replace('[O.]', '[O]')
            clean_smiles = clean_smiles.replace('[OO.]', 'OO')
            clean_smiles = clean_smiles.replace('O[O.]', 'O[O]')  # peroxy radical
            clean_smiles = clean_smiles.replace('(O.)', '([O])')
            clean_smiles = clean_smiles.replace('(OO.)', '(OO)')
            clean_smiles = clean_smiles.replace('[N+]([O-])', 'N(=O)')
            # Remove invalid patterns - any [X.] becomes [X]
            clean_smiles = re.sub(r'\[([A-Z][a-z]?)\.\]', r'[\1]', clean_smiles)

        # Format compatible with frontend: {id, label, smiles, raw_formula, code, title}
        nodes.append({
            "id": code,
            "label": species.formula[:20] if species.formula else code,
            "smiles": clean_smiles,
            "raw_formula": species.formula,
            "code": code,
            "title": f"{code}: {species.formula}",
            "molecular_weight": species.molecular_weight,
            "functional_groups": species.functional_groups,
            "is_radical": species.is_radical
        })

        if code in tree.pathways:
            for pathway in tree.pathways[code]:
                if pathway.branching_ratio < min_branching:
                    continue
                for product in pathway.products:
                    # First add the product node (if possible)
                    add_node(product, depth + 1)

                    # Only add edge if BOTH endpoints are in our node set
                    if product in node_ids:
                        product_species = tree.species_info.get(product, Species(code=product))
                        edges.append({
                            "from": code,
                            "to": product,
                            "yield": pathway.branching_ratio,
                            "reaction": f"{species.formula} + {pathway.oxidant} -> {product_species.formula}",
                            "label": f"{pathway.branching_ratio:.2f}" if pathway.branching_ratio != 1.0 else "",
                            "oxidant": pathway.oxidant
                        })

    add_node(tree.root, 0)

    # Output format matching what frontend expects: {"nodes": [...], "edges": [...]}
    result = {"nodes": nodes, "edges": edges}

    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)

    logger.info(f"Generated Cytoscape JSON: {len(nodes)} nodes, {len(edges)} edges")

    return result


def generate_kpp_format(tree: MechanismTree, output_path: str):
    """Generate KPP (Kinetics PreProcessor) format mechanism."""
    lines = []

    lines.append("// GECKO-A Generated Mechanism for KPP")
    lines.append(f"// Root VOC: {tree.root}")
    lines.append("")

    # Species list
    species_set = set()
    for reaction in tree.all_reactions:
        for s, _ in reaction.reactants:
            species_set.add(s)
        for s, _ in reaction.products:
            species_set.add(s)

    lines.append("#DEFVAR")
    for s in sorted(species_set):
        lines.append(f"  {s} = IGNORE;")
    lines.append("#ENDDEVAL")
    lines.append("")

    # Reactions
    lines.append("#EQUATIONS")
    for i, reaction in enumerate(tree.all_reactions):
        lhs = " + ".join([f"{c:.2f} {s}" if c != 1.0 else s for s, c in reaction.reactants])
        rhs = " + ".join([f"{c:.2f} {s}" if c != 1.0 else s for s, c in reaction.products])
        rate = reaction.get_rate_expression("kpp")
        lines.append(f"{{ {i + 1} }} {lhs} = {rhs} : {rate} ;")
    lines.append("#ENDEQUATIONS")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def generate_mcm_format(tree: MechanismTree, output_path: str):
    """Generate MCM (Master Chemical Mechanism) format."""
    lines = []

    lines.append("* GECKO-A Generated Mechanism for MCM")
    lines.append(f"* Root VOC: {tree.root}")
    lines.append("*")
    lines.append("* Generic Rate Coefficients")
    lines.append("*")

    for i, reaction in enumerate(tree.all_reactions):
        lhs = " + ".join([s for s, _ in reaction.reactants])
        rhs = " + ".join([s for s, _ in reaction.products])
        rate = reaction.get_rate_expression("mcm")
        lines.append(f"% {rate} : {lhs} = {rhs} ;")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def generate_facsimile_format(tree: MechanismTree, output_path: str):
    """Generate FACSIMILE format mechanism."""
    lines = []

    lines.append("* GECKO-A Generated Mechanism for FACSIMILE")
    lines.append(f"* Root VOC: {tree.root}")
    lines.append("*")
    lines.append("COMPILE INSTANT ;")
    lines.append("")
    lines.append("PARAMETER ;")
    lines.append("  TEMP = 298.0 ;")
    lines.append("")

    # Species
    species_set = set()
    for reaction in tree.all_reactions:
        for s, _ in reaction.reactants:
            species_set.add(s)
        for s, _ in reaction.products:
            species_set.add(s)

    lines.append("VARIABLE ;")
    for s in sorted(species_set):
        lines.append(f"  {s} ;")
    lines.append("")

    # Reactions
    for i, reaction in enumerate(tree.all_reactions):
        lhs = " + ".join([s for s, _ in reaction.reactants])
        rhs = " + ".join([s for s, _ in reaction.products])
        lines.append(f"% K{i + 1} : {lhs} = {rhs} ;")

    lines.append("")
    lines.append("COMPILE ;")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def _smiles_to_mol_image(smiles: str, size: Tuple[int, int] = (200, 150),
                          compound_name: Optional[str] = None,
                          formula: Optional[str] = None) -> Optional[Any]:
    """
    Convert SMILES to a molecule image with optional name and formula labels.

    Args:
        smiles: SMILES string
        size: Tuple of (width, height) for the image
        compound_name: Optional compound name to display below structure
        formula: Optional molecular formula to display below name

    Returns:
        numpy array of the image for matplotlib
    """
    if not HAS_RDKIT or not smiles:
        return None

    try:
        # Clean up SMILES - fix common radical notations for RDKit compatibility
        clean_smiles = smiles
        # Replace radical dot notations with standard SMILES
        # These patterns appear in GECKO atmospheric chemistry species
        clean_smiles = clean_smiles.replace('[O.]', '[O]')
        clean_smiles = clean_smiles.replace('[OO.]', 'O[O]')
        clean_smiles = clean_smiles.replace('O[O.]', 'O[O]')  # peroxy radical
        clean_smiles = clean_smiles.replace('(O.)', '([O])')
        clean_smiles = clean_smiles.replace('(OO.)', '(O[O])')
        # Also handle patterns at end of SMILES
        if clean_smiles.endswith('[O.]'):
            clean_smiles = clean_smiles[:-4] + '[O]'
        # Use regex to catch any remaining [X.] patterns
        clean_smiles = re.sub(r'\[([A-Z][a-z]?)\.\]', r'[\1]', clean_smiles)

        mol = Chem.MolFromSmiles(clean_smiles)
        if mol is None:
            # Try without stereochemistry
            clean_smiles = clean_smiles.replace('/', '').replace('\\', '')
            mol = Chem.MolFromSmiles(clean_smiles)

        if mol is None:
            return None

        # Get molecular formula if not provided
        if not formula:
            try:
                formula = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                pass

        # Generate 2D coordinates with high quality
        if hasattr(rdDepictor, 'SetPreferCoordGen'):
            rdDepictor.SetPreferCoordGen(True)

        # Determine if we should prepare molecule (RDKit 2020.03+)
        if hasattr(rdMolDraw2D, 'PrepareMolForDrawing'):
            rdMolDraw2D.PrepareMolForDrawing(mol)

        AllChem.Compute2DCoords(mol)

        # Calculate dimensions - reserve space for labels if needed
        label_height = 30 if (compound_name or formula) else 0
        mol_height = size[1] - label_height

        # Draw molecule
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], mol_height)
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        opts.addAtomIndices = False
        opts.bondLineWidth = 2.5
        opts.multipleBondOffset = 0.15
        opts.padding = 0.08
        opts.minFontSize = 12
        opts.annotationFontScale = 0.8

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Convert to PIL Image
        from io import BytesIO
        from PIL import Image, ImageDraw, ImageFont
        img_data = drawer.GetDrawingText()
        mol_img = Image.open(BytesIO(img_data))

        # If we have labels, create a combined image
        if compound_name or formula:
            # Create new image with space for labels
            combined_img = Image.new('RGB', (size[0], size[1]), 'white')
            combined_img.paste(mol_img, (0, 0))

            # Add labels
            draw = ImageDraw.Draw(combined_img)

            # Try to get a good font
            try:
                font_paths = [
                    '/System/Library/Fonts/Helvetica.ttc',
                    '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
                    'C:/Windows/Fonts/arial.ttf',
                ]
                font = None
                for fp in font_paths:
                    if os.path.exists(fp):
                        font = ImageFont.truetype(fp, 9)
                        break
                if font is None:
                    font = ImageFont.load_default()
            except Exception:
                font = ImageFont.load_default()

            y_pos = mol_height + 2

            # Draw compound name
            if compound_name:
                display_name = compound_name.replace('_', ' ')
                if len(display_name) > 18:
                    display_name = display_name[:16] + '..'
                bbox = draw.textbbox((0, 0), display_name, font=font)
                text_width = bbox[2] - bbox[0]
                x_pos = (size[0] - text_width) // 2
                draw.text((x_pos, y_pos), display_name, fill='#1a237e', font=font)
                y_pos += 12

            # Draw formula
            if formula:
                bbox = draw.textbbox((0, 0), formula, font=font)
                text_width = bbox[2] - bbox[0]
                x_pos = (size[0] - text_width) // 2
                draw.text((x_pos, y_pos), formula, fill='#424242', font=font)

            return np.array(combined_img)
        else:
            return np.array(mol_img)

    except Exception as e:
        logger.debug(f"Failed to render SMILES {smiles}: {e}")
        return None


def generate_static_diagram(tree: MechanismTree, output_path: str,
                           max_depth: int = 3, max_nodes: int = 25):
    """Generate static PNG/SVG mechanism diagram with RDKit molecular structures."""
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for static diagram")
        return

    # Build nodes to display
    nodes_by_depth: Dict[int, List[str]] = defaultdict(list)
    edges: List[Tuple[str, str, str, float]] = []  # (parent, child, oxidant, ratio)
    visited = set()

    def collect_nodes(code: str, depth: int):
        if code in visited or len(visited) >= max_nodes:
            return
        if depth > max_depth:
            return

        visited.add(code)
        nodes_by_depth[depth].append(code)

        if code in tree.pathways:
            for pathway in tree.pathways[code][:6]:  # Max 6 pathways per node
                for product in pathway.products[:2]:  # Max 2 products per pathway
                    edges.append((code, product, pathway.oxidant, pathway.branching_ratio))
                    collect_nodes(product, depth + 1)

    collect_nodes(tree.root, 0)

    # Calculate figure dimensions
    max_width = max(len(codes) for codes in nodes_by_depth.values()) if nodes_by_depth else 1
    n_depths = len(nodes_by_depth)

    # Fixed node size in inches
    node_w_inches = 1.5
    node_h_inches = 1.2
    h_spacing = 0.5  # horizontal spacing between nodes
    v_spacing = 1.0  # vertical spacing between levels

    fig_width = max(12, max_width * (node_w_inches + h_spacing))
    fig_height = max(10, n_depths * (node_h_inches + v_spacing) + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Position nodes (in inches from bottom-left)
    positions = {}
    for depth, codes in nodes_by_depth.items():
        y = fig_height - 1.5 - depth * (node_h_inches + v_spacing)
        n = len(codes)
        total_width = n * node_w_inches + (n - 1) * h_spacing
        start_x = (fig_width - total_width) / 2

        for i, code in enumerate(codes):
            x = start_x + i * (node_w_inches + h_spacing) + node_w_inches / 2
            positions[code] = (x, y)

    # Pre-render molecule images with names and formulas
    mol_images = {}
    if HAS_RDKIT:
        for code in positions:
            species = tree.species_info.get(code, Species(code=code))
            if species.smiles:
                # Get a human-readable display name using the mapping
                display_name = get_compound_name(code)
                # If we still have just the code, try to make it more readable
                if display_name == code:
                    # Truncate long codes
                    display_name = code[:12] + '..' if len(code) > 14 else code

                # Get molecular formula
                mol_formula = species.formula if species.formula else None
                # If formula is too long (GECKO format), try to extract molecular formula
                if mol_formula and len(mol_formula) > 15:
                    # Try to compute from SMILES instead
                    try:
                        test_mol = Chem.MolFromSmiles(species.smiles)
                        if test_mol:
                            mol_formula = rdMolDescriptors.CalcMolFormula(test_mol)
                    except Exception:
                        mol_formula = mol_formula[:15] + '..'

                img = _smiles_to_mol_image(
                    species.smiles,
                    (180, 160),  # Slightly taller to accommodate labels
                    compound_name=display_name,
                    formula=mol_formula
                )
                if img is not None:
                    mol_images[code] = img

    # Draw edges first
    drawn_edges = set()
    for parent, child, oxidant, ratio in edges:
        if parent in positions and child in positions:
            edge_key = (parent, child)
            if edge_key in drawn_edges:
                continue
            drawn_edges.add(edge_key)

            px, py = positions[parent]
            cx, cy = positions[child]

            # Arrow from bottom of parent to top of child
            ax.annotate(
                "",
                xy=(cx, cy + node_h_inches/2 + 0.05),
                xytext=(px, py - node_h_inches/2 - 0.05),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color='#333333',
                    lw=1.2,
                    connectionstyle="arc3,rad=0.15" if abs(px - cx) > 0.5 else "arc3,rad=0"
                )
            )

            # Oxidant label
            if oxidant:
                mid_x = (px + cx) / 2
                mid_y = (py - node_h_inches/2 + cy + node_h_inches/2) / 2
                label = f"+ {oxidant}"
                if ratio < 1.0 and ratio > 0:
                    label = f"{ratio:.2f}"
                ax.text(mid_x + 0.15, mid_y, label, fontsize=7, color='#555',
                       ha='left', va='center')

    # Draw nodes
    for code, (x, y) in positions.items():
        species = tree.species_info.get(code, Species(code=code))

        # Node rectangle position (bottom-left corner)
        rect_x = x - node_w_inches / 2
        rect_y = y - node_h_inches / 2

        if code in mol_images:
            # Draw molecule image using imshow with extent
            img = mol_images[code]
            ax.imshow(img, extent=[rect_x, rect_x + node_w_inches, rect_y, rect_y + node_h_inches],
                     aspect='auto', zorder=2)
            # Draw border
            rect = plt.Rectangle((rect_x, rect_y), node_w_inches, node_h_inches,
                                 fill=False, edgecolor='black', linewidth=1, zorder=3)
            ax.add_patch(rect)
        else:
            # Fallback: draw colored box with text
            color = '#ADD8E6' if species.is_radical else '#FFFACD'
            rect = FancyBboxPatch(
                (rect_x, rect_y), node_w_inches, node_h_inches,
                boxstyle="round,pad=0.02",
                facecolor=color,
                edgecolor='black',
                linewidth=1,
                zorder=2
            )
            ax.add_patch(rect)

            # Label
            label = species.formula[:18] if species.formula else code
            ax.text(x, y, label, ha='center', va='center', fontsize=7,
                   family='monospace', zorder=3)

    # Title - use human-readable name if available
    root_name = get_compound_name(tree.root)
    ax.set_title(f"Oxidation Mechanism: {root_name}", fontsize=14, fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')

    # Also save SVG
    try:
        plt.savefig(output_path.replace('.png', '.svg'), format='svg', bbox_inches='tight')
    except Exception as e:
        logger.warning(f"Failed to save SVG: {e}")

    plt.close()


# ==============================================================================
# Main Interface
# ==============================================================================

def generate_all_outputs(
    output_dir: str,
    voc_name: str,
    generate_static: bool = True,
    generate_kpp: bool = True,
    generate_mcm: bool = True,
    generate_facsimile: bool = True,
    max_depth: int = 3,
    min_branching: float = 0.05,
    max_nodes: int = 25,
    max_children: int = 4
) -> Dict[str, str]:
    """
    Generate all mechanism outputs.

    Returns dict of output paths.
    """
    results = {}

    # Parse mechanism
    parser = MechanismParser(output_dir, voc_name)
    tree = parser.parse()

    # Generate Cytoscape JSON (always)
    json_path = os.path.join(output_dir, "reaction_tree.json")
    cyto_data = generate_cytoscape_json(tree, json_path, max_depth, max_nodes, min_branching)
    results["reaction_tree"] = json_path

    # Save mechanism summary
    summary = {
        "root": tree.root,
        "root_formula": tree.root_formula,
        "total_species": len(tree.species_info),
        "total_reactions": len(tree.all_reactions),
        "pathways_count": sum(len(p) for p in tree.pathways.values())
    }
    summary_path = os.path.join(output_dir, "mechanism_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    results["summary"] = summary_path

    # Static diagram
    if generate_static and HAS_MATPLOTLIB:
        diagram_path = os.path.join(output_dir, "mechanism_diagram.png")
        generate_static_diagram(tree, diagram_path, max_depth, max_nodes)
        results["diagram_png"] = diagram_path
        results["diagram_svg"] = diagram_path.replace('.png', '.svg')

    # KPP format
    if generate_kpp:
        kpp_path = os.path.join(output_dir, f"{voc_name.lower()}.kpp")
        generate_kpp_format(tree, kpp_path)
        results["kpp"] = kpp_path

    # MCM format
    if generate_mcm:
        mcm_path = os.path.join(output_dir, f"{voc_name.lower()}.mcm")
        generate_mcm_format(tree, mcm_path)
        results["mcm"] = mcm_path

    # FACSIMILE format
    if generate_facsimile:
        fac_path = os.path.join(output_dir, f"{voc_name.lower()}.fac")
        generate_facsimile_format(tree, fac_path)
        results["facsimile"] = fac_path

    logger.info(f"Generated mechanism outputs: {list(results.keys())}")
    return results
