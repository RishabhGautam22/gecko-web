"""
GECKO-A Publication-Quality Pathway Visualizer

Generates hierarchical reaction pathway diagrams similar to GECKO-A publications,
featuring:
- Embedded molecular structure drawings (via RDKit)
- Hierarchical tree layout (via Graphviz)
- Branching ratio labels on edges
- Reagent annotations (+OH, +HO2/NO/RO2, etc.)
- Legend with variable definitions

Author: Deeksha Sharma
"""

import os
import re
import tempfile
import shutil
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Check dependencies
try:
    import graphviz
    HAS_GRAPHVIZ = True
except ImportError:
    HAS_GRAPHVIZ = False
    logger.warning("graphviz not available - pathway visualization disabled")

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor, Descriptors, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.warning("RDKit not available - will use text labels instead of structures")

# PIL for adding labels to images
try:
    from PIL import Image, ImageDraw, ImageFont
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    logger.warning("PIL not available - molecule labels may be limited")

# Import the unified SMILES conversion function and KNOWN_SPECIES database from reaction_tree
# This is the AUTHORITATIVE source for GECKO -> SMILES conversion
try:
    from gecko_web.reaction_tree import KNOWN_SPECIES, gecko_to_smiles as unified_gecko_to_smiles
    from gecko_web.reaction_tree import get_compound_name, GECKO_CODE_TO_NAME
except ImportError:
    KNOWN_SPECIES = {}
    unified_gecko_to_smiles = None
    GECKO_CODE_TO_NAME = {}
    def get_compound_name(code): return code


@dataclass
class SpeciesNode:
    """Represents a chemical species in the reaction pathway."""
    node_id: str
    smiles: str
    name: Optional[str] = None
    formula: Optional[str] = None
    is_product: bool = False
    is_radical: bool = False
    functional_groups: List[str] = field(default_factory=list)
    molecular_weight: Optional[float] = None


@dataclass
class ReactionEdge:
    """Represents a reaction step between species."""
    source: str
    target: str
    branching_ratio: Optional[float] = None
    reagents: List[str] = field(default_factory=list)
    oxidant: Optional[str] = None
    reaction_type: Optional[str] = None  # e.g., "addition", "abstraction", "isomerization"
    isomerization_percent: Optional[float] = None


class GECKOPathwayVisualizer:
    """
    Generate publication-quality reaction pathway diagrams
    compatible with GECKO-A atmospheric chemistry mechanisms.

    Produces diagrams similar to Figure 2 in GECKO-A publications,
    showing hierarchical oxidation pathways with molecular structures.
    """

    def __init__(self,
                 mol_size: Tuple[int, int] = (180, 140),
                 font_size: int = 11,
                 edge_font_size: int = 9,
                 rankdir: str = 'TB',
                 dpi: int = 150):
        """
        Initialize the pathway visualizer.

        Args:
            mol_size: Size of molecular structure images (width, height) in pixels
            font_size: Font size for node labels
            edge_font_size: Font size for edge labels
            rankdir: Graph direction ('TB' = top-bottom, 'LR' = left-right)
            dpi: Resolution for output images
        """
        if not HAS_GRAPHVIZ:
            raise ImportError("graphviz package is required for pathway visualization")

        self.nodes: Dict[str, SpeciesNode] = {}
        self.edges: List[ReactionEdge] = []
        self.mol_size = mol_size
        self.font_size = font_size
        self.edge_font_size = edge_font_size
        self.rankdir = rankdir
        self.dpi = dpi
        self.temp_dir = tempfile.mkdtemp(prefix='gecko_pathway_')
        self._image_cache: Dict[str, str] = {}

        # Graph for analysis
        if HAS_NETWORKX:
            self.graph = nx.DiGraph()
        else:
            self.graph = None

    def add_species(self,
                    node_id: str,
                    smiles: str,
                    name: Optional[str] = None,
                    formula: Optional[str] = None,
                    is_product: bool = False,
                    is_radical: bool = False,
                    functional_groups: Optional[List[str]] = None,
                    molecular_weight: Optional[float] = None) -> None:
        """
        Add a species node to the pathway.

        Args:
            node_id: Unique identifier for the species
            smiles: SMILES string for molecular structure
            name: Optional common name
            formula: Optional molecular formula
            is_product: Whether this is a terminal product
            is_radical: Whether this is a radical species
            functional_groups: List of functional groups present
            molecular_weight: Molecular weight in g/mol
        """
        node = SpeciesNode(
            node_id=node_id,
            smiles=smiles,
            name=name,
            formula=formula,
            is_product=is_product,
            is_radical=is_radical,
            functional_groups=functional_groups or [],
            molecular_weight=molecular_weight
        )
        self.nodes[node_id] = node

        if self.graph is not None:
            self.graph.add_node(node_id, **{
                'smiles': smiles,
                'name': name,
                'formula': formula,
                'is_product': is_product,
                'is_radical': is_radical
            })

    def add_reaction(self,
                     source: str,
                     target: str,
                     branching_ratio: Optional[float] = None,
                     reagents: Optional[List[str]] = None,
                     oxidant: Optional[str] = None,
                     reaction_type: Optional[str] = None,
                     isomerization_percent: Optional[float] = None) -> None:
        """
        Add a reaction edge between species.

        Args:
            source: Source species node_id
            target: Target species node_id
            branching_ratio: Branching ratio (0-1)
            reagents: List of reagents (e.g., ['OH'], ['HO2', 'NO', 'RO2'])
            oxidant: Primary oxidant for this step
            reaction_type: Type of reaction
            isomerization_percent: Percentage for isomerization reactions
        """
        edge = ReactionEdge(
            source=source,
            target=target,
            branching_ratio=branching_ratio,
            reagents=reagents or [],
            oxidant=oxidant,
            reaction_type=reaction_type,
            isomerization_percent=isomerization_percent
        )
        self.edges.append(edge)

        if self.graph is not None:
            self.graph.add_edge(source, target, **{
                'branching_ratio': branching_ratio,
                'reagents': reagents or [],
                'oxidant': oxidant
            })

    def _clean_smiles_for_rdkit(self, smiles: str, node_id: str = "") -> str:
        """
        Clean SMILES string for RDKit compatibility.
        Handles GECKO-specific notation.

        Args:
            smiles: Original SMILES string
            node_id: Optional node ID for KNOWN_SPECIES lookup
        """
        if not smiles and not node_id:
            return ""

        # PRESERVE OVERRIDES: If the incoming SMILES is explicitly set to Isopentane structure,
        # do NOT let the dictionary lookup revert it to n-pentane. This handles the fix
        # applied in reaction_tree.py based on job name.
        if smiles == 'CC(C)CC' and node_id and ('05000' in node_id or 'PENTANE' in node_id.upper()):
            return smiles

        # First, try to get a known good SMILES from the database
        if node_id and KNOWN_SPECIES:
            node_upper = node_id.upper().strip()
            for variant in [node_upper, node_upper.lstrip('G')]:
                if variant in KNOWN_SPECIES:
                    return KNOWN_SPECIES[variant]

        # If SMILES looks like invalid aromatic notation (common GECKO issue)
        # e.g., 'c()ccccc0C' is invalid - use unified conversion
        if smiles and ('c()' in smiles or 'c0' in smiles.lower() or 'c1H' in smiles):
            if unified_gecko_to_smiles and node_id:
                better_smiles = unified_gecko_to_smiles("", species_code=node_id)
                if better_smiles and 'c()' not in better_smiles and better_smiles != smiles:
                    logger.info(f"Using unified SMILES for {node_id}: {better_smiles}")
                    return better_smiles

        if not smiles:
            return ""

        clean = smiles

        # Remove radical dots that RDKit doesn't understand
        clean = re.sub(r'\[([A-Z][a-z]?)\.\]', r'[\1]', clean)
        clean = clean.replace('[O.]', '[O]')
        clean = clean.replace('[OO.]', 'OO')
        clean = clean.replace('(O.)', '([O])')
        clean = clean.replace('(OO.)', '(OO)')

        # Handle peroxy radicals
        clean = clean.replace('O[O.]', 'O[O]')
        clean = clean.replace('[N+]([O-])', 'N(=O)')

        # Fix nitrate groups - ONO2 should be ON(=O)=O but sometimes appears incorrectly
        clean = clean.replace('ON(O)=O', 'ON(=O)=O')
        clean = clean.replace('N(O)=O', 'N(=O)=O')

        # Fix cumulated double bonds (allene-like) that shouldn't exist
        # e.g., C=C=C is valid (allene) but CC(=C=C)=O is likely an error
        # These often come from incorrect GECKO->SMILES conversion
        # Replace =C=C with =CC for terminal methylenes
        if '=C=C)' in clean:
            clean = clean.replace('=C=C)', '=CC)')

        # Remove trailing dots
        clean = re.sub(r'\.$', '', clean)

        return clean

    def _generate_mol_image(self,
                            smiles: str,
                            node_id: str,
                            highlight_radical: bool = False,
                            compound_name: Optional[str] = None,
                            formula: Optional[str] = None) -> Optional[str]:
        """
        Generate molecular structure image from SMILES with name and formula labels.

        Args:
            smiles: SMILES string for the molecule
            node_id: Unique identifier for this node
            highlight_radical: Whether to highlight radical sites
            compound_name: Optional compound name to display below structure
            formula: Optional molecular formula to display below name

        Returns path to the generated PNG file.
        """
        if not HAS_RDKIT:
            return self._generate_text_placeholder(node_id)

        # Check cache - include name and formula in cache key
        cache_key = f"{smiles}_{node_id}_{compound_name}_{formula}"
        if cache_key in self._image_cache:
            return self._image_cache[cache_key]

        clean_smiles = self._clean_smiles_for_rdkit(smiles, node_id=node_id)

        # Early validation - if SMILES is empty or too short, use placeholder
        if not clean_smiles or len(clean_smiles) < 1:
            logger.warning(f"Empty SMILES for {node_id}, using placeholder")
            return self._generate_text_placeholder(node_id, node_id)

        try:
            mol = Chem.MolFromSmiles(clean_smiles)

            # If parsing fails, try without sanitization as a last resort
            if mol is None:
                mol = Chem.MolFromSmiles(clean_smiles, sanitize=False)
                if mol is not None:
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception as san_e:
                        # Sanitization failed - molecule is likely invalid
                        # Use placeholder instead of rendering broken structure
                        logger.warning(f"SMILES sanitization failed for {node_id}: {san_e}")
                        return self._generate_text_placeholder(node_id, node_id)

            if mol is None:
                logger.warning(f"Could not parse SMILES for {node_id}: {clean_smiles}")
                return self._generate_text_placeholder(node_id, node_id)

            # Additional validation - check if molecule has atoms
            if mol.GetNumAtoms() == 0:
                logger.warning(f"Empty molecule for {node_id}")
                return self._generate_text_placeholder(node_id, node_id)

            # Get molecular formula if not provided
            if not formula and HAS_RDKIT:
                try:
                    formula = rdMolDescriptors.CalcMolFormula(mol)
                except Exception:
                    pass

            # Compute 2D coordinates with high quality
            try:
                if hasattr(rdDepictor, 'SetPreferCoordGen'):
                    rdDepictor.SetPreferCoordGen(True)

                # Prepare molecule for drawing (canonicalize, add stereo, etc)
                if hasattr(rdMolDraw2D, 'PrepareMolForDrawing'):
                    rdMolDraw2D.PrepareMolForDrawing(mol)

                AllChem.Compute2DCoords(mol)
            except Exception as coord_e:
                logger.warning(f"Could not compute 2D coords for {node_id}: {coord_e}")
                return self._generate_text_placeholder(node_id, node_id)

            # Calculate image dimensions - add space for labels
            label_height = 40 if (compound_name or formula) else 0
            mol_width = self.mol_size[0]
            mol_height = self.mol_size[1] - label_height  # Reserve space for labels

            # Create high-quality drawing (Cairo/PNG)
            drawer = rdMolDraw2D.MolDraw2DCairo(mol_width, mol_height)

            # Drawing options for publication quality
            opts = drawer.drawOptions()
            opts.addStereoAnnotation = True
            opts.addAtomIndices = False
            opts.bondLineWidth = 2.5  # Thicker bonds for clarity
            opts.multipleBondOffset = 0.15
            opts.padding = 0.08  # Slightly more padding to accommodate labels
            opts.minFontSize = 12 # Ensure labels are readable
            opts.annotationFontScale = 0.8

            # Draw the molecule
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()

            # Save molecule image to temporary file
            mol_img_data = drawer.GetDrawingText()

            # If PIL is available and we have labels, add them below the structure
            if HAS_PIL and (compound_name or formula):
                img_path = self._add_labels_to_mol_image(
                    mol_img_data, node_id, compound_name, formula, mol_width, mol_height, label_height
                )
            else:
                # Just save the molecule image without labels
                img_path = os.path.join(self.temp_dir, f"{node_id}.png")
                with open(img_path, 'wb') as f:
                    f.write(mol_img_data)

            self._image_cache[cache_key] = img_path
            return img_path

        except Exception as e:
            logger.warning(f"Error generating structure for {node_id}: {e}")
            return self._generate_text_placeholder(node_id, node_id)

    def _add_labels_to_mol_image(self,
                                  mol_img_data: bytes,
                                  node_id: str,
                                  compound_name: Optional[str],
                                  formula: Optional[str],
                                  mol_width: int,
                                  mol_height: int,
                                  label_height: int) -> str:
        """
        Add compound name and formula labels below the molecule structure.

        Args:
            mol_img_data: PNG image data from RDKit
            node_id: Node identifier for filename
            compound_name: Compound name to display
            formula: Molecular formula to display
            mol_width: Width of molecule image
            mol_height: Height of molecule image
            label_height: Height reserved for labels

        Returns:
            Path to the labeled image file
        """
        import io

        # Load the molecule image
        mol_img = Image.open(io.BytesIO(mol_img_data))

        # Create a new image with space for labels
        total_height = mol_height + label_height
        combined_img = Image.new('RGB', (mol_width, total_height), 'white')

        # Paste molecule image at top
        combined_img.paste(mol_img, (0, 0))

        # Draw labels
        draw = ImageDraw.Draw(combined_img)

        # Try to use a nice font, fall back to default
        try:
            # Try common font paths
            font_paths = [
                '/System/Library/Fonts/Helvetica.ttc',  # macOS
                '/System/Library/Fonts/SFNSText.ttf',   # macOS
                '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',  # Linux
                '/usr/share/fonts/TTF/DejaVuSans.ttf',  # Arch Linux
                'C:/Windows/Fonts/arial.ttf',  # Windows
            ]
            font_name = None
            font_formula = None
            for fp in font_paths:
                if os.path.exists(fp):
                    font_name = ImageFont.truetype(fp, 11)
                    font_formula = ImageFont.truetype(fp, 10)
                    break
            if font_name is None:
                font_name = ImageFont.load_default()
                font_formula = ImageFont.load_default()
        except Exception:
            font_name = ImageFont.load_default()
            font_formula = ImageFont.load_default()

        # Position for labels (centered below molecule)
        y_pos = mol_height + 2

        # Draw compound name (bold-ish by drawing twice with slight offset)
        if compound_name:
            # Truncate long names
            display_name = compound_name.replace('_', ' ')
            if len(display_name) > 20:
                display_name = display_name[:18] + '...'

            # Get text bounding box for centering
            bbox = draw.textbbox((0, 0), display_name, font=font_name)
            text_width = bbox[2] - bbox[0]
            x_pos = (mol_width - text_width) // 2

            # Draw name in dark blue
            draw.text((x_pos, y_pos), display_name, fill='#1a237e', font=font_name)
            y_pos += 14

        # Draw molecular formula
        if formula:
            # Get text bounding box for centering
            bbox = draw.textbbox((0, 0), formula, font=font_formula)
            text_width = bbox[2] - bbox[0]
            x_pos = (mol_width - text_width) // 2

            # Draw formula in gray
            draw.text((x_pos, y_pos), formula, fill='#424242', font=font_formula)

        # Save combined image
        img_path = os.path.join(self.temp_dir, f"{node_id}_labeled.png")
        combined_img.save(img_path, 'PNG')

        return img_path

    def _generate_text_placeholder(self,
                                   node_id: str,
                                   formula: Optional[str] = None) -> str:
        """Generate a text-based placeholder when structure can't be drawn."""
        # Create a simple SVG placeholder
        display_text = formula[:20] if formula else node_id
        if len(display_text) > 20:
            display_text = display_text[:17] + "..."

        svg_content = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{self.mol_size[0]}" height="{self.mol_size[1]}">
    <rect width="100%" height="100%" fill="#e8f4f8" stroke="#4a90a4" stroke-width="2" rx="5"/>
    <text x="50%" y="45%" text-anchor="middle"
          font-family="Arial, sans-serif" font-size="11" font-weight="bold" fill="#333">
        {node_id}
    </text>
    <text x="50%" y="65%" text-anchor="middle"
          font-family="Arial, sans-serif" font-size="9" fill="#666">
        {display_text}
    </text>
</svg>'''

        svg_path = os.path.join(self.temp_dir, f"{node_id}_placeholder.svg")
        with open(svg_path, 'w') as f:
            f.write(svg_content)
        return svg_path

    def _format_edge_label(self, edge: ReactionEdge) -> str:
        """
        Format edge label with branching ratio and reagents.
        Follows GECKO-A publication style.
        """
        parts = []

        # Branching ratio
        if edge.branching_ratio is not None and edge.branching_ratio < 1.0:
            parts.append(f"{edge.branching_ratio:.3f}")

        # Isomerization percentage
        if edge.isomerization_percent is not None:
            parts.append(f"{edge.isomerization_percent:.0f}%")

        # Reagents - format like "+ OH" or "+ HO2/NO/RO2"
        if edge.reagents:
            if len(edge.reagents) == 1:
                parts.append(f"+ {edge.reagents[0]}")
            else:
                # Group related reagents
                reagent_str = "+ " + "/".join(edge.reagents)
                parts.append(reagent_str)
        elif edge.oxidant:
            parts.append(f"+ {edge.oxidant}")

        return "\\n".join(parts) if parts else ""

    def _calculate_edge_style(self, edge: ReactionEdge) -> Dict[str, str]:
        """Calculate edge styling based on branching ratio."""
        style = {
            'color': '#333333',
            'penwidth': '1.5',
            'arrowsize': '0.8'
        }

        if edge.branching_ratio is not None:
            # Thicker lines for major pathways
            width = max(1.0, min(3.0, edge.branching_ratio * 4))
            style['penwidth'] = f'{width:.1f}'

            # Color based on importance
            if edge.branching_ratio >= 0.3:
                style['color'] = '#000000'  # Black for major pathways
            elif edge.branching_ratio >= 0.1:
                style['color'] = '#333333'  # Dark gray
            else:
                style['color'] = '#666666'  # Light gray for minor pathways

        return style

    def generate_diagram(self,
                         output_path: str,
                         format: str = 'png',
                         title: Optional[str] = None,
                         legend: Optional[Dict[str, str]] = None,
                         show_legend: bool = True) -> str:
        """
        Generate the complete reaction pathway diagram.

        Args:
            output_path: Output file path (without extension)
            format: Output format ('png', 'svg', 'pdf')
            title: Optional diagram title
            legend: Optional legend dictionary {symbol: description}
            show_legend: Whether to show the legend

        Returns:
            Path to the generated diagram file
        """
        # Create Graphviz digraph
        dot = graphviz.Digraph(
            comment='GECKO-A Reaction Pathway',
            format=format,
            engine='dot'
        )

        # Graph-level attributes for publication quality
        dot.attr(
            rankdir=self.rankdir,
            splines='ortho',  # Use orthogonal edges for cleaner layout
            nodesep='0.8',
            ranksep='1.5',
            fontname='Arial',
            fontsize=str(self.font_size),
            bgcolor='white',
            dpi=str(self.dpi),
            overlap='false',
            concentrate='false',  # Do not merge edges - improves label readability
            margin='0.5'  # Add margin to prevent title/nodes from being clipped
        )

        # Node defaults
        dot.attr('node',
                 shape='box',
                 style='filled',
                 fillcolor='white',
                 fontname='Arial',
                 fontsize=str(self.font_size))

        # Edge defaults
        dot.attr('edge',
                 fontname='Arial',
                 fontsize=str(self.edge_font_size),
                 arrowhead='normal')

        # Add title
        if title:
            dot.attr(label=f'<<B>{title}</B>>',
                     labelloc='t',
                     fontsize=str(self.font_size + 4))

        # Generate molecule images and add nodes
        for node_id, node in self.nodes.items():
            # Pass compound name and formula for labeling below structure
            img_path = self._generate_mol_image(
                node.smiles,
                node_id,
                node.is_radical,
                compound_name=node.name,
                formula=node.formula
            )

            if img_path and os.path.exists(img_path):
                # Use image as node
                node_attrs = {
                    'label': '',
                    'shape': 'none',
                    'image': img_path,
                    'imagescale': 'true',
                    'fixedsize': 'true',
                    'width': str(self.mol_size[0] / 72),  # Convert px to inches
                    'height': str(self.mol_size[1] / 72)
                }
            else:
                # Fallback to text label
                label = node.name or node_id
                if node.formula:
                    label = f"{label}\\n{node.formula[:30]}"
                node_attrs = {
                    'label': label,
                    'shape': 'box',
                    'style': 'filled,rounded',
                    'fillcolor': '#e8f4f8' if not node.is_radical else '#ffe4e4'
                }

            # Mark products
            if node.is_product:
                node_attrs['xlabel'] = '<<I>products</I>>'

            dot.node(node_id, **node_attrs)

        # Add edges with labels and styling
        for edge in self.edges:
            label = self._format_edge_label(edge)
            style = self._calculate_edge_style(edge)

            edge_attrs = {
                'label': f" {label} ",  # Add spaces for padding
                'color': style['color'],
                'penwidth': style['penwidth'],
                'arrowsize': style['arrowsize'],
                'fontcolor': '#333333' # Dark gray for contrast
            }

            dot.edge(edge.source, edge.target, **edge_attrs)

        # Add legend as a subgraph
        if show_legend and legend:
            with dot.subgraph(name='cluster_legend') as legend_graph:
                legend_graph.attr(label='Legend', fontsize=str(self.font_size))
                legend_graph.attr(style='dashed', color='gray')

                legend_text = "\\l".join([f"{k} = {v}" for k, v in legend.items()]) + "\\l"
                legend_graph.node('legend_content',
                                  label=legend_text,
                                  shape='box',
                                  style='filled',
                                  fillcolor='#f9f9f9',
                                  fontsize=str(self.font_size - 1))

        # Render the diagram
        try:
            output_file = dot.render(output_path, cleanup=True)
            logger.info(f"Generated pathway diagram: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Error rendering diagram: {e}")
            raise

    def from_reaction_tree(self, tree_data: Dict) -> None:
        """
        Build the pathway from a reaction tree dictionary.

        Args:
            tree_data: Dictionary with 'nodes' and 'edges' lists
        """
        # Add nodes
        for node in tree_data.get('nodes', []):
            node_code = node.get('id') or node.get('code')
            # Get human-readable name for the compound
            display_name = get_compound_name(node_code)
            self.add_species(
                node_id=node_code,
                smiles=node.get('smiles', ''),
                name=display_name,
                formula=node.get('raw_formula') or node.get('label'),
                is_radical=node.get('is_radical', False),
                functional_groups=node.get('functional_groups', []),
                molecular_weight=node.get('molecular_weight')
            )

        # Add edges
        for edge in tree_data.get('edges', []):
            reagents = []
            if edge.get('oxidant'):
                reagents.append(edge['oxidant'])

            self.add_reaction(
                source=edge.get('from') or edge.get('source'),
                target=edge.get('to') or edge.get('target'),
                branching_ratio=edge.get('yield') or edge.get('branching_ratio'),
                reagents=reagents,
                oxidant=edge.get('oxidant')
            )

    def from_gecko_mechanism(self,
                             mechanism_file: str,
                             dictionary_file: str,
                             root_species: str,
                             max_depth: int = 3,
                             min_branching: float = 0.05) -> None:
        """
        Parse GECKO-A mechanism files and build the pathway.

        Args:
            mechanism_file: Path to reactions.txt or .mec file
            dictionary_file: Path to dictionary.out file
            root_species: Starting species code (e.g., 'APINEN', 'ISOPRE')
            max_depth: Maximum reaction depth to include
            min_branching: Minimum branching ratio to include
        """
        # Parse dictionary for species info
        species_info = {}
        if os.path.exists(dictionary_file):
            with open(dictionary_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        code = parts[0]
                        formula = parts[1]
                        mw = float(parts[2]) if len(parts) > 2 else None
                        species_info[code] = {
                            'formula': formula,
                            'molecular_weight': mw
                        }

        # Parse mechanism file for reactions
        reactions = []
        with open(mechanism_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('!'):
                    continue

                # Parse reaction: REACTANT + OXIDANT -> PRODUCTS
                if '->' in line or '=' in line:
                    reactions.append(line)

        # Build adjacency list
        adj_list = {}
        for rxn in reactions:
            # Simple parsing - this should be enhanced for complex mechanisms
            try:
                if '->' in rxn:
                    lhs, rhs = rxn.split('->')
                elif '=' in rxn:
                    lhs, rhs = rxn.split('=')
                else:
                    continue

                reactants = [r.strip() for r in lhs.split('+')]
                products = [p.strip() for p in rhs.split('+')]

                # Find organic reactant (not OH, HO2, NO, etc.)
                inorganic = {'OH', 'HO2', 'NO', 'NO2', 'NO3', 'O3', 'O2', 'H2O', 'CO', 'CO2'}
                organic_reactant = None
                oxidant = None

                for r in reactants:
                    r_clean = r.split()[0] if ' ' in r else r
                    if r_clean in inorganic:
                        oxidant = r_clean
                    else:
                        organic_reactant = r_clean

                if organic_reactant:
                    if organic_reactant not in adj_list:
                        adj_list[organic_reactant] = []

                    for p in products:
                        p_clean = p.split()[0] if ' ' in p else p
                        if p_clean not in inorganic:
                            adj_list[organic_reactant].append({
                                'to': p_clean,
                                'oxidant': oxidant,
                                'yield': 1.0  # Default yield
                            })
            except Exception as e:
                logger.debug(f"Could not parse reaction: {rxn} - {e}")
                continue

        # BFS to build tree from root
        visited = set()
        queue = [(root_species, 0)]

        while queue:
            current, depth = queue.pop(0)

            if current in visited or depth > max_depth:
                continue

            visited.add(current)

            # Add species node
            info = species_info.get(current, {})
            self.add_species(
                node_id=current,
                smiles=self._gecko_to_smiles(info.get('formula', current), species_code=current),
                formula=info.get('formula'),
                molecular_weight=info.get('molecular_weight')
            )

            # Add reactions
            if current in adj_list:
                for child_info in adj_list[current]:
                    child = child_info['to']
                    yield_val = child_info.get('yield', 1.0)

                    if yield_val < min_branching:
                        continue

                    reagents = [child_info['oxidant']] if child_info.get('oxidant') else []

                    self.add_reaction(
                        source=current,
                        target=child,
                        branching_ratio=yield_val,
                        reagents=reagents,
                        oxidant=child_info.get('oxidant')
                    )

                    if child not in visited:
                        queue.append((child, depth + 1))

    def _gecko_to_smiles(self, gecko_formula: str, species_code: str = "") -> str:
        """
        Convert GECKO formula notation to SMILES.
        Delegates to the unified gecko_to_smiles function in reaction_tree.py.

        Args:
            gecko_formula: GECKO-style formula notation
            species_code: Optional species code for database lookup

        Returns:
            Valid SMILES string
        """
        # Use the unified conversion function from reaction_tree
        if unified_gecko_to_smiles:
            return unified_gecko_to_smiles(gecko_formula, species_code)

        # Fallback if import failed (shouldn't happen in normal operation)
        if not gecko_formula:
            return ""
        return gecko_formula

    def cleanup(self) -> None:
        """Clean up temporary files."""
        if os.path.exists(self.temp_dir):
            try:
                shutil.rmtree(self.temp_dir)
            except Exception as e:
                logger.warning(f"Could not cleanup temp dir: {e}")

    def __del__(self):
        """Destructor to ensure cleanup."""
        self.cleanup()


def generate_pathway_diagram(tree_data: Dict,
                             output_path: str,
                             title: Optional[str] = None,
                             format: str = 'png',
                             mol_size: Tuple[int, int] = (180, 140)) -> str:
    """
    Convenience function to generate a pathway diagram from tree data.

    Args:
        tree_data: Dictionary with 'nodes' and 'edges' lists
        output_path: Output file path (without extension)
        title: Optional diagram title
        format: Output format ('png', 'svg', 'pdf')
        mol_size: Size of molecular structure images

    Returns:
        Path to the generated diagram file
    """
    viz = GECKOPathwayVisualizer(mol_size=mol_size)
    viz.from_reaction_tree(tree_data)

    # Default legend for atmospheric chemistry
    legend = {
        'X': '-OOH, -ONO2, -OH',
        'RO2': 'Peroxy radicals'
    }

    output_file = viz.generate_diagram(
        output_path=output_path,
        format=format,
        title=title,
        legend=legend
    )

    viz.cleanup()
    return output_file


# Test function
if __name__ == '__main__':
    # Test with a simple example
    test_tree = {
        'nodes': [
            {'id': 'APINEN', 'smiles': 'CC1=CCC2CC1C2(C)C', 'code': 'APINEN', 'formula': 'C10H16'},
            {'id': 'PROD1', 'smiles': 'CC1(C)C2CCC(O)C(C)C2C1', 'code': 'PROD1', 'is_radical': True},
            {'id': 'PROD2', 'smiles': 'CC1(C)C2CC(O)CC(C)C2C1', 'code': 'PROD2'}
        ],
        'edges': [
            {'from': 'APINEN', 'to': 'PROD1', 'yield': 0.439, 'oxidant': 'OH'},
            {'from': 'APINEN', 'to': 'PROD2', 'yield': 0.220, 'oxidant': 'OH'}
        ]
    }

    output = generate_pathway_diagram(
        test_tree,
        'test_pathway',
        title='Test: Alpha-Pinene + OH Oxidation',
        format='png'
    )
    print(f"Generated: {output}")
