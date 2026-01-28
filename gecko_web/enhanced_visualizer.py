"""
GECKO-A Enhanced Pathway Visualizer with High-Quality RDKit Rendering

This module provides publication-quality reaction pathway diagrams with:
- Server-side RDKit rendering for molecular structures
- Reaction type color coding (OH addition, H-abstraction, isomerization, etc.)
- Branching ratio annotations with variable line weights
- Variable legend support (X = -OOH, -ONO2, -OH)
- Multiple layout algorithms (hierarchical, radial, force-directed)
- Stereochemistry handling for E/Z isomers
- Interactive SVG output with hover information

References:
- GECKO-A oxidation mechanism diagrams (Aumont et al.)
- MCM diagram conventions
- IUPAC atmospheric chemistry nomenclature

Author: GECKO-A Development Team
"""

import os
import re
import logging
import tempfile
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from enum import Enum
import base64

logger = logging.getLogger(__name__)


# ==============================================================================
# Dependency Checks
# ==============================================================================

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdDepictor
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.warning("RDKit not available - molecular rendering disabled")

try:
    import graphviz
    HAS_GRAPHVIZ = True
except ImportError:
    HAS_GRAPHVIZ = False
    logger.warning("graphviz not available - diagram generation disabled")



# ==============================================================================
# Reaction Type Classification and Colors
# ==============================================================================

class ReactionType(Enum):
    """Classification of atmospheric oxidation reaction types."""
    OH_ADDITION = "oh_addition"
    OH_ABSTRACTION = "oh_abstraction"
    O3_OZONOLYSIS = "o3_ozonolysis"
    NO3_ADDITION = "no3_addition"
    PHOTOLYSIS = "photolysis"
    ISOMERIZATION = "isomerization"
    DECOMPOSITION = "decomposition"
    RO2_NO = "ro2_no"
    RO2_HO2 = "ro2_ho2"
    RO2_RO2 = "ro2_ro2"
    UNKNOWN = "unknown"


# Color scheme for reaction types (publication-quality colors)
REACTION_TYPE_COLORS = {
    ReactionType.OH_ADDITION: "#0072B2",      # Blue - OH addition
    ReactionType.OH_ABSTRACTION: "#009E73",   # Green - H-abstraction
    ReactionType.O3_OZONOLYSIS: "#D55E00",    # Orange - Ozonolysis
    ReactionType.NO3_ADDITION: "#CC79A7",     # Pink - NO3 reactions
    ReactionType.PHOTOLYSIS: "#F0E442",       # Yellow - Photolysis
    ReactionType.ISOMERIZATION: "#56B4E9",    # Light blue - Isomerization
    ReactionType.DECOMPOSITION: "#E69F00",    # Gold - Decomposition
    ReactionType.RO2_NO: "#999999",           # Gray - RO2+NO
    ReactionType.RO2_HO2: "#AAAAAA",          # Light gray - RO2+HO2
    ReactionType.RO2_RO2: "#BBBBBB",          # Lighter gray - RO2+RO2
    ReactionType.UNKNOWN: "#666666",          # Dark gray - Unknown
}

# Functional group markers for variable legend
FUNCTIONAL_GROUPS = {
    'OOH': {'pattern': 'OO', 'name': 'Hydroperoxide', 'marker': 'X'},
    'ONO2': {'pattern': 'O[N+](=O)[O-]', 'name': 'Nitrate', 'marker': 'X'},
    'OH': {'pattern': 'O', 'name': 'Hydroxyl', 'marker': 'X'},
    'OO': {'pattern': 'O[O]', 'name': 'Peroxy radical', 'marker': '●'},
    'CHO': {'pattern': 'C=O', 'name': 'Aldehyde', 'marker': '△'},
    'CO': {'pattern': 'C(=O)C', 'name': 'Ketone', 'marker': '◇'},
    'COOH': {'pattern': 'C(=O)O', 'name': 'Carboxylic acid', 'marker': '□'},
    'epoxy': {'pattern': 'C1OC1', 'name': 'Epoxide', 'marker': '○'},
}


# ==============================================================================
# Layout Algorithms
# ==============================================================================

class LayoutAlgorithm(Enum):
    """Available layout algorithms for diagrams."""
    HIERARCHICAL = "hierarchical"  # Top-down tree layout (default)
    RADIAL = "radial"              # Circular/radial layout
    FORCE_DIRECTED = "force"       # Physics-based layout
    KAMADA_KAWAI = "kamada"        # Energy minimization
    CUSTOM = "custom"              # User-defined positions


# ==============================================================================
# Data Classes
# ==============================================================================

@dataclass
class MoleculeRenderOptions:
    """Options for molecular structure rendering."""
    width: int = 200
    height: int = 150
    bond_line_width: float = 2.0
    atom_label_font_size: int = 14
    show_atom_numbers: bool = False
    add_stereo_annotation: bool = True
    highlight_stereo: bool = True
    highlight_radicals: bool = True
    radical_color: Tuple[float, float, float] = (1.0, 0.0, 0.0)  # Red
    background_color: Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)
    use_svg: bool = False
    wedge_bonds: bool = True


@dataclass
class DiagramStyle:
    """Style options for the complete diagram."""
    # Graph layout
    rankdir: str = 'TB'  # TB, LR, BT, RL
    ranksep: float = 1.2  # Vertical separation
    nodesep: float = 0.8  # Horizontal separation

    # Fonts
    title_font_size: int = 16
    node_font_size: int = 11
    edge_font_size: int = 9
    legend_font_size: int = 10

    # Colors
    background_color: str = '#FFFFFF'
    node_border_color: str = '#333333'
    edge_color: str = '#666666'

    # Edge styling
    edge_min_width: float = 1.0
    edge_max_width: float = 4.0
    use_reaction_colors: bool = True
    show_branching_ratios: bool = True
    color_by_reaction_type: bool = True

    # Output
    dpi: int = 300
    format: str = 'png'  # png, svg, pdf


@dataclass
class SpeciesData:
    """Complete data for a chemical species node."""
    code: str
    smiles: str
    formula: Optional[str] = None
    name: Optional[str] = None
    molecular_weight: Optional[float] = None
    n_carbon: int = 0
    n_oxygen: int = 0
    n_nitrogen: int = 0
    is_radical: bool = False
    is_product: bool = False
    functional_groups: List[str] = field(default_factory=list)
    position: Optional[Tuple[float, float]] = None
    depth: int = 0


@dataclass
class ReactionData:
    """Complete data for a reaction edge."""
    source: str
    target: str
    branching_ratio: float = 1.0
    reaction_type: ReactionType = ReactionType.UNKNOWN
    oxidant: Optional[str] = None
    reagents: List[str] = field(default_factory=list)
    isomerization_fraction: Optional[float] = None
    rate_constant: Optional[float] = None


# ==============================================================================
# Molecular Renderer
# ==============================================================================

class MolecularRenderer:
    """
    High-quality molecular structure renderer using RDKit.

    Features:
    - Publication-quality 2D depictions
    - Stereochemistry display (wedge bonds, E/Z annotations)
    - Radical site highlighting
    - Functional group highlighting
    """

    def __init__(self, options: Optional[MoleculeRenderOptions] = None):
        if not HAS_RDKIT:
            raise ImportError("RDKit is required for molecular rendering")

        self.options = options or MoleculeRenderOptions()
        self._cache: Dict[str, bytes] = {}

    def render_smiles(self, smiles: str,
                      highlight_atoms: Optional[List[int]] = None,
                      highlight_bonds: Optional[List[int]] = None) -> Optional[bytes]:
        """
        Render a SMILES string to an image.

        Args:
            smiles: SMILES string
            highlight_atoms: Optional list of atom indices to highlight
            highlight_bonds: Optional list of bond indices to highlight

        Returns:
            PNG or SVG image data as bytes, or None on failure
        """
        cache_key = f"{smiles}_{highlight_atoms}_{highlight_bonds}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        try:
            mol = self._parse_smiles(smiles)
            if mol is None:
                return None

            # Compute 2D coordinates
            rdDepictor.Compute2DCoords(mol)

            # Create drawer
            if self.options.use_svg:
                drawer = rdMolDraw2D.MolDraw2DSVG(
                    self.options.width, self.options.height
                )
            else:
                drawer = rdMolDraw2D.MolDraw2DCairo(
                    self.options.width, self.options.height
                )

            # Configure drawing options
            self._configure_draw_options(drawer)

            # Prepare highlighting
            highlight_atom_map = {}
            highlight_bond_map = {}

            if highlight_atoms:
                highlight_atom_colors = {
                    i: self.options.radical_color for i in highlight_atoms
                }
                highlight_atom_map = highlight_atom_colors

            # Draw molecule
            if highlight_atoms or highlight_bonds:
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms or [],
                    highlightAtomColors=highlight_atom_map or {},
                    highlightBonds=highlight_bonds or [],
                    highlightBondColors=highlight_bond_map or {}
                )
            else:
                drawer.DrawMolecule(mol)

            drawer.FinishDrawing()

            img_data = drawer.GetDrawingText()
            if isinstance(img_data, str):
                img_data = img_data.encode('utf-8')

            self._cache[cache_key] = img_data
            return img_data

        except Exception as e:
            logger.warning(f"Failed to render SMILES '{smiles}': {e}")
            return None

    def _parse_smiles(self, smiles: str) -> Optional[Any]:
        """Parse SMILES with cleanup for GECKO-specific notation."""
        # Clean up radical notation
        clean = smiles
        clean = re.sub(r'\[([A-Z][a-z]?)\.\]', r'[\1]', clean)
        clean = clean.replace('[O.]', '[O]')
        clean = clean.replace('(O.)', '([O])')
        clean = clean.replace('O[O.]', 'O[O]')

        try:
            mol = Chem.MolFromSmiles(clean)
            if mol is not None:
                return mol

            # Try without sanitization
            mol = Chem.MolFromSmiles(clean, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except:
                    pass
                return mol

        except Exception:
            pass

        return None

    def _configure_draw_options(self, drawer):
        """Configure RDKit drawing options for publication quality."""
        opts = drawer.drawOptions()
        opts.bondLineWidth = self.options.bond_line_width
        opts.addStereoAnnotation = self.options.add_stereo_annotation
        opts.addAtomIndices = self.options.show_atom_numbers
        opts.padding = 0.1
        opts.multipleBondOffset = 0.15

        # Set background
        if len(self.options.background_color) >= 3:
            opts.setBackgroundColour(self.options.background_color)

    def detect_radical_sites(self, smiles: str) -> List[int]:
        """Detect radical sites in a molecule for highlighting."""
        radical_atoms = []

        try:
            mol = self._parse_smiles(smiles)
            if mol is None:
                return []

            for atom in mol.GetAtoms():
                # Check for unpaired electrons
                if atom.GetNumRadicalElectrons() > 0:
                    radical_atoms.append(atom.GetIdx())
                # Also check for [O] pattern (common in GECKO)
                if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == 0:
                    if atom.GetDegree() == 1:
                        radical_atoms.append(atom.GetIdx())

        except Exception:
            pass

        return radical_atoms

    def detect_stereocenters(self, smiles: str) -> Dict[str, Any]:
        """Detect stereocenters and E/Z bonds in a molecule."""
        stereo_info = {
            'chiral_centers': [],
            'ez_bonds': [],
            'has_stereochemistry': False
        }

        try:
            mol = self._parse_smiles(smiles)
            if mol is None:
                return stereo_info

            # Detect chiral centers
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            stereo_info['chiral_centers'] = chiral_centers

            # Detect E/Z bonds
            for bond in mol.GetBonds():
                stereo = bond.GetStereo()
                if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
                    stereo_info['ez_bonds'].append({
                        'bond_idx': bond.GetIdx(),
                        'stereo': 'E' if stereo == Chem.BondStereo.STEREOE else 'Z',
                        'atoms': (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    })

            stereo_info['has_stereochemistry'] = bool(
                chiral_centers or stereo_info['ez_bonds']
            )

        except Exception:
            pass

        return stereo_info

    def render_to_base64(self, smiles: str) -> Optional[str]:
        """Render SMILES to base64-encoded image string."""
        img_data = self.render_smiles(smiles)
        if img_data:
            return base64.b64encode(img_data).decode('ascii')
        return None


# ==============================================================================
# Reaction Type Classifier
# ==============================================================================

class ReactionClassifier:
    """Classifies atmospheric oxidation reactions by type."""

    @staticmethod
    def classify_reaction(source_smiles: str,
                         target_smiles: str,
                         oxidant: Optional[str] = None,
                         reagents: Optional[List[str]] = None) -> ReactionType:
        """
        Classify a reaction based on reactants, products, and oxidant.

        Args:
            source_smiles: SMILES of reactant
            target_smiles: SMILES of product
            oxidant: Primary oxidant (OH, O3, NO3, etc.)
            reagents: Additional reagents

        Returns:
            ReactionType enum value
        """
        oxidant_upper = (oxidant or "").upper()
        reagents_upper = [r.upper() for r in (reagents or [])]

        # Check for specific oxidants
        if oxidant_upper == "OH" or "OH" in reagents_upper:
            # Distinguish addition vs abstraction
            if ReactionClassifier._is_addition(source_smiles, target_smiles):
                return ReactionType.OH_ADDITION
            return ReactionType.OH_ABSTRACTION

        if oxidant_upper == "O3" or "O3" in reagents_upper:
            return ReactionType.O3_OZONOLYSIS

        if oxidant_upper == "NO3" or "NO3" in reagents_upper:
            return ReactionType.NO3_ADDITION

        if "HV" in reagents_upper or "PHOT" in oxidant_upper:
            return ReactionType.PHOTOLYSIS

        # Check for RO2 reactions
        if "NO" in reagents_upper and "RO2" not in oxidant_upper:
            return ReactionType.RO2_NO

        if "HO2" in reagents_upper:
            return ReactionType.RO2_HO2

        if "RO2" in reagents_upper:
            return ReactionType.RO2_RO2

        # Check for isomerization (same formula, different structure)
        if ReactionClassifier._is_isomerization(source_smiles, target_smiles):
            return ReactionType.ISOMERIZATION

        return ReactionType.UNKNOWN

    @staticmethod
    def _is_addition(source_smiles: str, target_smiles: str) -> bool:
        """Check if reaction is an addition (adds atoms to skeleton)."""
        if not HAS_RDKIT:
            return False

        try:
            source_mol = Chem.MolFromSmiles(source_smiles)
            target_mol = Chem.MolFromSmiles(target_smiles)

            if source_mol is None or target_mol is None:
                return False

            # Addition typically increases heavy atom count
            source_heavy = source_mol.GetNumHeavyAtoms()
            target_heavy = target_mol.GetNumHeavyAtoms()

            return target_heavy > source_heavy

        except Exception:
            return False

    @staticmethod
    def _is_isomerization(source_smiles: str, target_smiles: str) -> bool:
        """Check if reaction is an isomerization."""
        if not HAS_RDKIT:
            return False

        try:
            source_mol = Chem.MolFromSmiles(source_smiles)
            target_mol = Chem.MolFromSmiles(target_smiles)

            if source_mol is None or target_mol is None:
                return False

            # Isomerization: same formula, different structure
            source_formula = Chem.rdMolDescriptors.CalcMolFormula(source_mol)
            target_formula = Chem.rdMolDescriptors.CalcMolFormula(target_mol)

            return source_formula == target_formula

        except Exception:
            return False


# ==============================================================================
# Enhanced Pathway Visualizer
# ==============================================================================

class EnhancedPathwayVisualizer:
    """
    Publication-quality pathway diagram generator with enhanced features.

    Key improvements:
    - High-resolution RDKit molecular rendering
    - Reaction type color coding
    - Branching ratio visualization
    - Variable legend support
    - Multiple layout algorithms
    - Stereochemistry handling
    """

    def __init__(self,
                 style: Optional[DiagramStyle] = None,
                 render_options: Optional[MoleculeRenderOptions] = None):
        """
        Initialize the enhanced visualizer.

        Args:
            style: Diagram styling options
            render_options: Molecular rendering options
        """
        self.style = style or DiagramStyle()
        self.render_options = render_options or MoleculeRenderOptions()

        self.species: Dict[str, SpeciesData] = {}
        self.reactions: List[ReactionData] = []

        if HAS_RDKIT:
            self.renderer = MolecularRenderer(self.render_options)
        else:
            self.renderer = None

        self.classifier = ReactionClassifier()
        self._temp_dir = tempfile.mkdtemp(prefix='gecko_viz_')
        self._image_files: Dict[str, str] = {}

    def add_species(self, species: SpeciesData) -> None:
        """Add a species to the diagram."""
        self.species[species.code] = species

    def add_reaction(self, reaction: ReactionData) -> None:
        """Add a reaction to the diagram."""
        # Auto-classify if not specified
        if reaction.reaction_type == ReactionType.UNKNOWN:
            source = self.species.get(reaction.source)
            target = self.species.get(reaction.target)

            if source and target:
                reaction.reaction_type = self.classifier.classify_reaction(
                    source.smiles,
                    target.smiles,
                    reaction.oxidant,
                    reaction.reagents
                )

        self.reactions.append(reaction)

    def from_tree_data(self, tree_data: Dict) -> None:
        """
        Load pathway from reaction tree JSON data.

        Args:
            tree_data: Dictionary with 'nodes' and 'edges' keys
        """
        # Load species
        for node in tree_data.get('nodes', []):
            code = node.get('id') or node.get('code', 'UNKNOWN')
            smiles = node.get('smiles', '')

            # Detect functional groups
            func_groups = node.get('functional_groups', [])
            if not func_groups:
                 # Fallback to internal detection if missing
                 func_groups = []
                 for fg_name, fg_info in FUNCTIONAL_GROUPS.items():
                     if fg_info['pattern'] in smiles:
                         func_groups.append(fg_name)

            # Calculate properties if RDKit available
            mw = node.get('molecular_weight')
            if mw is None and HAS_RDKIT and smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        mw = Descriptors.MolWt(mol)
                except:
                    pass

            species = SpeciesData(
                code=code,
                smiles=smiles,
                formula=node.get('raw_formula') or node.get('label'),
                name=node.get('name') or node.get('code'),
                molecular_weight=mw,
                is_radical=node.get('is_radical',
                                   'O[O]' in smiles or '[O]' in smiles),
                is_product=node.get('is_product', False),
                functional_groups=func_groups
            )
            self.add_species(species)

        # Load reactions
        for edge in tree_data.get('edges', []):
            source = edge.get('from') or edge.get('source')
            target = edge.get('to') or edge.get('target')

            reagents = []
            oxidant = edge.get('oxidant')
            if oxidant:
                reagents.append(oxidant)

            reaction = ReactionData(
                source=source,
                target=target,
                branching_ratio=edge.get('yield') or edge.get('branching_ratio', 1.0),
                oxidant=oxidant,
                reagents=reagents
            )
            self.add_reaction(reaction)

    def generate(self,
                 output_path: str,
                 title: Optional[str] = None,
                 show_legend: bool = True,
                 layout: LayoutAlgorithm = LayoutAlgorithm.HIERARCHICAL) -> str:
        """
        Generate the pathway diagram.

        Args:
            output_path: Output file path (without extension)
            title: Optional diagram title
            show_legend: Whether to show the variable legend
            layout: Layout algorithm to use

        Returns:
            Path to generated diagram file
        """
        if not HAS_GRAPHVIZ:
            raise ImportError("graphviz is required for diagram generation")

        # Generate molecular images
        self._generate_molecular_images()

        # Create graphviz diagram
        dot = graphviz.Digraph(
            comment='GECKO-A Reaction Pathway',
            format=self.style.format,
            engine=self._get_layout_engine(layout)
        )

        # Graph attributes
        dot.attr(
            rankdir=self.style.rankdir,
            splines='ortho' if layout == LayoutAlgorithm.HIERARCHICAL else 'polyline',
            nodesep=str(self.style.nodesep),
            ranksep=str(self.style.ranksep),
            bgcolor=self.style.background_color,
            dpi=str(self.style.dpi),
            fontname='Arial',
            fontsize=str(self.style.node_font_size)
        )

        # Add title
        if title:
            dot.attr(
                label=f'<<B>{title}</B>>',
                labelloc='t',
                fontsize=str(self.style.title_font_size)
            )

        # Add nodes
        for code, species in self.species.items():
            self._add_node(dot, species)

        # Add edges
        for reaction in self.reactions:
            self._add_edge(dot, reaction)

        # Add legend
        if show_legend:
            self._add_legend(dot)

        # Render
        try:
            output_file = dot.render(output_path, cleanup=True)
            logger.info(f"Generated diagram: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Failed to render diagram: {e}")
            raise

    def generate_diagram(self,
                         output_path: str,
                         title: Optional[str] = None,
                         layout: str = "hierarchical",
                         show_branching_ratios: bool = True,
                         color_by_reaction_type: bool = True,
                         format: str = "png") -> str:
        """
        Generate pathway diagram with extended options (API-compatible method).

        Args:
            output_path: Output file path (without extension)
            title: Optional diagram title
            layout: Layout algorithm ('hierarchical', 'radial', 'force')
            show_branching_ratios: Whether to show branching ratios on edges
            color_by_reaction_type: Whether to color edges by reaction type
            format: Output format ('png', 'svg', 'pdf')

        Returns:
            Path to generated diagram file
        """
        # Map string layout to enum
        layout_map = {
            'hierarchical': LayoutAlgorithm.HIERARCHICAL,
            'radial': LayoutAlgorithm.RADIAL,
            'force': LayoutAlgorithm.FORCE_DIRECTED,
            'force-directed': LayoutAlgorithm.FORCE_DIRECTED,
        }
        layout_enum = layout_map.get(layout.lower(), LayoutAlgorithm.HIERARCHICAL)

        # Update style based on options
        self.style.format = format
        self.style.show_branching_ratios = show_branching_ratios
        self.style.color_by_reaction_type = color_by_reaction_type

        return self.generate(
            output_path=output_path,
            title=title,
            show_legend=True,
            layout=layout_enum
        )

    def _get_layout_engine(self, layout: LayoutAlgorithm) -> str:
        """Get Graphviz engine for layout algorithm."""
        engine_map = {
            LayoutAlgorithm.HIERARCHICAL: 'dot',
            LayoutAlgorithm.RADIAL: 'twopi',
            LayoutAlgorithm.FORCE_DIRECTED: 'neato',
            LayoutAlgorithm.KAMADA_KAWAI: 'fdp',
            LayoutAlgorithm.CUSTOM: 'neato',
        }
        return engine_map.get(layout, 'dot')

    def _generate_molecular_images(self) -> None:
        """Generate molecular structure images for all species."""
        for code, species in self.species.items():
            if not species.smiles:
                continue

            if self.renderer:
                # Detect radical sites for highlighting
                radical_sites = []
                if species.is_radical:
                    radical_sites = self.renderer.detect_radical_sites(species.smiles)

                img_data = self.renderer.render_smiles(
                    species.smiles,
                    highlight_atoms=radical_sites if radical_sites else None
                )

                if img_data:
                    img_path = os.path.join(self._temp_dir, f"{code}.png")
                    with open(img_path, 'wb') as f:
                        f.write(img_data)
                    self._image_files[code] = img_path

    def _add_node(self, dot: "graphviz.Digraph", species: SpeciesData) -> None:
        """Add a species node to the diagram."""
        node_attrs = {}

        # Check for image
        if species.code in self._image_files:
            img_path = self._image_files[species.code]
            node_attrs.update({
                'label': '',
                'shape': 'none',
                'image': img_path,
                'imagescale': 'true',
                'fixedsize': 'true',
                'width': str(self.render_options.width / 72),
                'height': str(self.render_options.height / 72)
            })
        else:
            # Text fallback
            label = species.name or species.code
            if species.formula:
                label += f"\\n{species.formula[:30]}"

            node_attrs.update({
                'label': label,
                'shape': 'box',
                'style': 'filled,rounded',
                'fillcolor': '#ffe4e4' if species.is_radical else '#e8f4f8',
                'fontsize': str(self.style.node_font_size)
            })

        # Add tooltip with additional info
        tooltip_parts = [species.code]
        if species.name:
            tooltip_parts.append(species.name)
        if species.molecular_weight:
            tooltip_parts.append(f"MW: {species.molecular_weight:.1f}")
        if species.functional_groups:
            tooltip_parts.append(f"Groups: {', '.join(species.functional_groups)}")

        node_attrs['tooltip'] = ' | '.join(tooltip_parts)

        dot.node(species.code, **node_attrs)

    def _add_edge(self, dot: "graphviz.Digraph", reaction: ReactionData) -> None:
        """Add a reaction edge to the diagram."""
        # Build edge label
        label_parts = []

        if reaction.branching_ratio < 1.0:
            label_parts.append(f"{reaction.branching_ratio:.3f}")

        if reaction.isomerization_fraction is not None:
            label_parts.append(f"{reaction.isomerization_fraction:.0f}%")

        if reaction.reagents:
            if len(reaction.reagents) == 1:
                label_parts.append(f"+ {reaction.reagents[0]}")
            else:
                label_parts.append(f"+ {'/'.join(reaction.reagents)}")
        elif reaction.oxidant:
            label_parts.append(f"+ {reaction.oxidant}")

        label = "\\n".join(label_parts)

        # Calculate edge styling
        edge_attrs = {
            'label': label,
            'fontsize': str(self.style.edge_font_size),
            'arrowhead': 'normal',
            'arrowsize': '0.8'
        }

        # Color by reaction type
        if self.style.use_reaction_colors:
            edge_attrs['color'] = REACTION_TYPE_COLORS.get(
                reaction.reaction_type, '#666666'
            )

        # Line width by branching ratio
        width_range = self.style.edge_max_width - self.style.edge_min_width
        width = self.style.edge_min_width + (reaction.branching_ratio * width_range)
        edge_attrs['penwidth'] = str(min(width, self.style.edge_max_width))

        dot.edge(reaction.source, reaction.target, **edge_attrs)

    def _add_legend(self, dot: "graphviz.Digraph") -> None:
        """Add variable and reaction type legend."""
        with dot.subgraph(name='cluster_legend') as legend:
            legend.attr(
                label='Legend',
                style='dashed',
                color='gray',
                fontsize=str(self.style.legend_font_size)
            )

            # Variable legend
            variable_text = "X = -OOH, -ONO2, -OH\\l"
            variable_text += "● = Peroxy radical\\l"
            variable_text += "△ = Aldehyde\\l"

            legend.node('legend_variables',
                       label=variable_text,
                       shape='box',
                       style='filled',
                       fillcolor='#f9f9f9',
                       fontsize=str(self.style.legend_font_size - 1))

            # Reaction type legend (if using colors)
            if self.style.use_reaction_colors:
                rxn_legend_text = ""
                for rxn_type, color in list(REACTION_TYPE_COLORS.items())[:6]:
                    name = rxn_type.value.replace('_', ' ').title()
                    rxn_legend_text += f"<FONT COLOR='{color}'>━━</FONT> {name}\\l"

                legend.node('legend_reactions',
                           label=f'<{rxn_legend_text}>',
                           shape='box',
                           style='filled',
                           fillcolor='#f9f9f9',
                           fontsize=str(self.style.legend_font_size - 1))

    def cleanup(self) -> None:
        """Clean up temporary files."""
        if os.path.exists(self._temp_dir):
            try:
                shutil.rmtree(self._temp_dir)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp dir: {e}")

    def __del__(self):
        self.cleanup()


# ==============================================================================
# Convenience Functions
# ==============================================================================

def generate_enhanced_pathway(
    tree_data: Dict,
    output_path: str,
    title: Optional[str] = None,
    format: str = 'png',
    mol_size: Tuple[int, int] = (200, 150),
    use_reaction_colors: bool = True,
    layout: str = 'hierarchical'
) -> str:
    """
    Generate an enhanced pathway diagram from tree data.

    Args:
        tree_data: Dictionary with 'nodes' and 'edges' keys
        output_path: Output file path (without extension)
        title: Optional diagram title
        format: Output format ('png', 'svg', 'pdf')
        mol_size: Size of molecular structures (width, height)
        use_reaction_colors: Color edges by reaction type
        layout: Layout algorithm ('hierarchical', 'radial', 'force')

    Returns:
        Path to generated diagram file
    """
    style = DiagramStyle(
        format=format,
        use_reaction_colors=use_reaction_colors
    )

    render_opts = MoleculeRenderOptions(
        width=mol_size[0],
        height=mol_size[1]
    )

    layout_enum = {
        'hierarchical': LayoutAlgorithm.HIERARCHICAL,
        'radial': LayoutAlgorithm.RADIAL,
        'force': LayoutAlgorithm.FORCE_DIRECTED,
    }.get(layout, LayoutAlgorithm.HIERARCHICAL)

    viz = EnhancedPathwayVisualizer(style=style, render_options=render_opts)
    viz.from_tree_data(tree_data)

    output = viz.generate(
        output_path=output_path,
        title=title,
        layout=layout_enum
    )

    viz.cleanup()
    return output


# ==============================================================================
# Module Test
# ==============================================================================

if __name__ == '__main__':
    # Test with sample data
    test_data = {
        'nodes': [
            {'id': 'APINEN', 'smiles': 'CC1=CCC2CC1C2(C)C', 'code': 'APINEN',
             'raw_formula': 'C10H16', 'name': 'α-Pinene'},
            {'id': 'PROD1', 'smiles': 'CC1(C)C2CCC(O)C(C)C2C1', 'code': 'PROD1',
             'is_radical': False},
            {'id': 'PROD2', 'smiles': 'CC1(C)C2CC([O])CC(C)C2C1', 'code': 'PROD2',
             'is_radical': True}
        ],
        'edges': [
            {'from': 'APINEN', 'to': 'PROD1', 'yield': 0.439, 'oxidant': 'OH'},
            {'from': 'APINEN', 'to': 'PROD2', 'yield': 0.220, 'oxidant': 'OH'}
        ]
    }

    output = generate_enhanced_pathway(
        test_data,
        'test_enhanced_pathway',
        title='Test: α-Pinene + OH Oxidation',
        format='png',
        use_reaction_colors=True
    )
    print(f"Generated: {output}")
