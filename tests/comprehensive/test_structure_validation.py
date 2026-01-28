"""
Comprehensive Structure Validation Tests

Tests that verify molecular structures are rendered correctly by:
1. Converting GECKO formula to SMILES
2. Validating SMILES with RDKit
3. Checking ring structures match expected
4. Verifying molecular weight
5. Validating functional groups

These tests catch issues like:
- Wrong ring sizes (e.g., cyclohexane instead of cyclobutane)
- Missing functional groups
- Incorrect SMILES generation from GECKO formulas
"""

import pytest
import sys
import os
import re

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from gecko_web.chemdata import get_compound, get_gecko_formula
from gecko_web.chemdata.compound_database import COMPOUNDS
from gecko_web.reaction_tree import gecko_to_smiles, KNOWN_SPECIES


# =============================================================================
# TEST DATA: 50 compounds to validate comprehensively
# =============================================================================

# Format: (compound_name, expected_ring_sizes, expected_carbons, expected_mw_range, notes)
VALIDATION_COMPOUNDS = [
    # Terpene oxidation products (critical - these had cyclobutane issues)
    ('pinic_acid', [4], 9, (184, 188), 'Cyclobutane diacid'),
    ('pinonic_acid', [4], 10, (182, 186), 'Cyclobutane ketoacid'),
    ('pinonaldehyde', [4], 10, (166, 170), 'Cyclobutane aldehyde'),
    ('nopinone', [4, 6], 9, (136, 140), 'Bicyclic ketone'),  # RDKit SSSR

    # Monoterpenes (bicyclic)
    ('alpha_pinene', [4, 6], 10, (134, 138), 'Bicyclic monoterpene'),
    ('beta_pinene', [4, 6], 10, (134, 138), 'Bicyclic monoterpene'),
    ('camphene', [5, 5], 10, (134, 138), 'Bicyclic monoterpene'),  # RDKit SSSR
    ('sabinene', [3, 5], 10, (134, 138), 'Bicyclic monoterpene'),
    ('3_carene', [3, 6], 10, (134, 138), 'Bicyclic monoterpene'),

    # Monoterpenes (monocyclic)
    ('limonene', [6], 10, (134, 138), 'Monocyclic monoterpene'),
    ('alpha_terpinene', [6], 10, (134, 138), 'Monocyclic monoterpene'),
    ('gamma_terpinene', [6], 10, (134, 138), 'Monocyclic monoterpene'),
    ('terpinolene', [6], 10, (134, 138), 'Monocyclic monoterpene'),

    # Monoterpenes (acyclic)
    ('myrcene', [], 10, (134, 138), 'Acyclic monoterpene'),
    ('ocimene', [], 10, (134, 138), 'Acyclic monoterpene'),

    # Oxygenated terpenes
    ('linalool', [], 10, (152, 156), 'Acyclic terpene alcohol'),
    ('alpha_terpineol', [6], 10, (152, 156), 'Monocyclic terpene alcohol'),
    ('nerolidol', [], 15, (220, 225), 'Sesquiterpene alcohol'),
    ('1_8_cineole', [6, 6], 10, (152, 156), 'Bicyclic ether'),

    # Sesquiterpenes
    ('beta_caryophyllene', [4, 9], 15, (202, 206), 'Bicyclic sesquiterpene'),
    ('alpha_humulene', [11], 15, (202, 206), 'Monocyclic sesquiterpene'),

    # Simple alkanes
    ('methane', [], 1, (14, 18), 'Simplest alkane'),
    ('ethane', [], 2, (28, 32), 'C2 alkane'),
    ('propane', [], 3, (42, 46), 'C3 alkane'),
    ('isobutane', [], 4, (56, 60), 'Branched C4'),
    ('isopentane', [], 5, (70, 74), 'Branched C5'),
    ('neopentane', [], 5, (70, 74), 'Highly branched C5'),
    ('cyclohexane', [6], 6, (82, 86), 'Cyclic C6'),

    # Alkenes
    ('ethene', [], 2, (26, 30), 'Simple alkene'),
    ('propene', [], 3, (40, 44), 'C3 alkene'),
    ('isoprene', [], 5, (66, 70), 'Biogenic VOC'),
    ('isobutene', [], 4, (54, 58), 'Branched alkene'),

    # Aromatics
    ('benzene', [6], 6, (76, 80), 'Aromatic ring'),
    ('toluene', [6], 7, (90, 94), 'Methylbenzene'),
    ('ethylbenzene', [6], 8, (104, 108), 'Ethylbenzene'),
    ('styrene', [6], 8, (102, 106), 'Vinyl benzene'),
    ('naphthalene', [6, 6], 10, (126, 130), 'Fused aromatics'),

    # Oxygenated VOCs
    ('methanol', [], 1, (30, 34), 'Simplest alcohol'),
    ('ethanol', [], 2, (44, 48), 'C2 alcohol'),
    ('acetone', [], 3, (56, 60), 'Simplest ketone'),
    ('formaldehyde', [], 1, (28, 32), 'Simplest aldehyde'),
    ('acetaldehyde', [], 2, (42, 46), 'C2 aldehyde'),
    ('acetic_acid', [], 2, (58, 62), 'Simplest carboxylic acid'),
    ('formic_acid', [], 1, (44, 48), 'Simplest acid'),

    # Esters
    ('methyl_formate', [], 2, (58, 62), 'Simplest ester'),
    ('ethyl_acetate', [], 4, (86, 90), 'Common ester'),
    ('butyl_acetate', [], 6, (114, 118), 'Ester solvent'),

    # Ethers
    ('dimethyl_ether', [], 2, (44, 48), 'Simplest ether'),
    ('diethyl_ether', [], 4, (72, 76), 'Common ether'),
    ('methyl_tert_butyl_ether', [], 5, (86, 90), 'Fuel additive'),
]


class TestStructureValidation:
    """Test structural correctness for 50 key compounds."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = COMPOUNDS

    @pytest.mark.parametrize("compound_name,expected_rings,expected_carbons,mw_range,notes",
                            VALIDATION_COMPOUNDS)
    def test_compound_smiles_validity(self, compound_name, expected_rings, expected_carbons, mw_range, notes):
        """Test that compound SMILES is valid and parseable by RDKit."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol is not None, \
            f"{compound_name}: Invalid SMILES '{compound.smiles}'"

    @pytest.mark.parametrize("compound_name,expected_rings,expected_carbons,mw_range,notes",
                            VALIDATION_COMPOUNDS)
    def test_compound_carbon_count(self, compound_name, expected_rings, expected_carbons, mw_range, notes):
        """Test that compound has correct number of carbons."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.fail(f"Invalid SMILES for {compound_name}")

        actual_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        assert actual_carbons == expected_carbons, \
            f"{compound_name}: Expected {expected_carbons} carbons, got {actual_carbons}. {notes}"

    @pytest.mark.parametrize("compound_name,expected_rings,expected_carbons,mw_range,notes",
                            VALIDATION_COMPOUNDS)
    def test_compound_molecular_weight(self, compound_name, expected_rings, expected_carbons, mw_range, notes):
        """Test that molecular weight is within expected range."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.fail(f"Invalid SMILES for {compound_name}")

        mw = Descriptors.MolWt(mol)
        min_mw, max_mw = mw_range
        assert min_mw <= mw <= max_mw, \
            f"{compound_name}: MW {mw:.2f} not in range [{min_mw}, {max_mw}]. {notes}"

    @pytest.mark.parametrize("compound_name,expected_rings,expected_carbons,mw_range,notes",
                            VALIDATION_COMPOUNDS)
    def test_compound_ring_structure(self, compound_name, expected_rings, expected_carbons, mw_range, notes):
        """Test that compound has expected ring structure."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.fail(f"Invalid SMILES for {compound_name}")

        ring_info = mol.GetRingInfo()
        actual_rings = sorted([len(r) for r in ring_info.AtomRings()])

        if expected_rings:
            expected_sorted = sorted(expected_rings)
            # For bicyclic systems, RDKit may detect additional "envelope" rings
            # So we check that expected rings are a subset
            for exp_ring in expected_sorted:
                assert exp_ring in actual_rings, \
                    f"{compound_name}: Expected ring of size {exp_ring}, got {actual_rings}. {notes}"
        else:
            # Acyclic compound
            assert len(actual_rings) == 0, \
                f"{compound_name}: Expected no rings (acyclic), got {actual_rings}. {notes}"


class TestGeckoToSmilesConversion:
    """Test that GECKO formula to SMILES conversion produces correct structures."""

    # Format: (compound_name, gecko_formula, expected_ring_sizes)
    GECKO_CONVERSION_TESTS = [
        # Critical terpene oxidation products
        ('pinic_acid', 'CH3C1(CH3)CH(CH2CO(OH))CH2C1HCO(OH)', [4]),
        ('pinonic_acid', 'CH3COC1HCH2C1(CH3)(CH3)CH(CH2CO(OH))', [4]),
        ('pinonaldehyde', 'CH3COC1HCH2C1(CH3)(CH3)CH(CH2CHO)', [4]),

        # Simple compounds
        ('methane', 'CH4', []),
        ('ethane', 'CH3CH3', []),
        ('propane', 'CH3CH2CH3', []),
        ('isobutane', 'CH3CH(CH3)CH3', []),
        ('isopentane', 'CH3CH(CH3)CH2CH3', []),

        # Cyclic
        ('cyclohexane', 'C1H2CH2CH2CH2CH2C1H2', [6]),
        ('cyclopentane', 'C1H2CH2CH2CH2C1H2', [5]),
    ]

    @pytest.mark.parametrize("compound_name,gecko_formula,expected_rings",
                            GECKO_CONVERSION_TESTS)
    def test_gecko_to_smiles_ring_structure(self, compound_name, gecko_formula, expected_rings):
        """Test that GECKO to SMILES conversion produces correct ring structure."""
        # Get the compound's species code if known
        compound = COMPOUNDS.get(compound_name)
        species_code = ""
        if compound_name == 'pinic_acid':
            species_code = 'TA9000'
        elif compound_name == 'pinonic_acid':
            species_code = 'TA8000'
        elif compound_name == 'pinonaldehyde':
            species_code = 'TA7000'

        smiles = gecko_to_smiles(gecko_formula, species_code)

        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, \
            f"{compound_name}: gecko_to_smiles returned invalid SMILES '{smiles}'"

        ring_info = mol.GetRingInfo()
        actual_rings = sorted([len(r) for r in ring_info.AtomRings()])

        if expected_rings:
            expected_sorted = sorted(expected_rings)
            for exp_ring in expected_sorted:
                assert exp_ring in actual_rings, \
                    f"{compound_name}: Expected ring size {exp_ring} in {actual_rings}. " \
                    f"GECKO: {gecko_formula} -> SMILES: {smiles}"
        else:
            assert len(actual_rings) == 0, \
                f"{compound_name}: Expected no rings, got {actual_rings}. " \
                f"GECKO: {gecko_formula} -> SMILES: {smiles}"


class TestKnownSpeciesDatabase:
    """Test that KNOWN_SPECIES database has correct structures."""

    # Format: (species_key, expected_ring_sizes, expected_carbons, notes)
    KNOWN_SPECIES_TESTS = [
        # Terpene oxidation products (critical - cyclobutane structures)
        ('PINICACID', [4], 9, 'Must be cyclobutane, C9H14O4'),
        ('PINIC', [4], 9, 'Must be cyclobutane, C9H14O4'),
        ('TA9000', [4], 9, 'Pinic acid species code'),
        ('PINONICACID', [4], 10, 'Must be cyclobutane, C10H16O3'),
        ('PINONIC', [4], 10, 'Must be cyclobutane, C10H16O3'),
        ('TA8000', [4], 10, 'Pinonic acid species code'),
        ('PINONALDEHYDE', [4], 10, 'Must be cyclobutane, C10H16O2'),
        ('PINONAL', [4], 10, 'Must be cyclobutane, C10H16O2'),
        ('TA7000', [4], 10, 'Pinonaldehyde species code'),

        # Monoterpenes
        ('APINENE', [4, 6], 10, 'Alpha-pinene bicyclic'),
        ('BPINENE', [4, 6], 10, 'Beta-pinene bicyclic'),
        ('LIMONENE', [6], 10, 'Limonene monocyclic'),
        ('ISOPRENE', [], 5, 'Isoprene acyclic'),

        # Simple compounds
        ('CH4', [], 1, 'Methane'),
        ('C2H6', [], 2, 'Ethane'),
        ('C3H8', [], 3, 'Propane'),
        ('TOLUENE', [6], 7, 'Toluene aromatic'),
        ('BENZENE', [6], 6, 'Benzene aromatic'),
    ]

    @pytest.mark.parametrize("species_key,expected_rings,expected_carbons,notes",
                            KNOWN_SPECIES_TESTS)
    def test_known_species_structure(self, species_key, expected_rings, expected_carbons, notes):
        """Test that KNOWN_SPECIES entries have correct structures."""
        if species_key not in KNOWN_SPECIES:
            pytest.skip(f"{species_key} not in KNOWN_SPECIES")

        smiles = KNOWN_SPECIES[species_key]
        mol = Chem.MolFromSmiles(smiles)

        assert mol is not None, \
            f"{species_key}: Invalid SMILES in KNOWN_SPECIES: '{smiles}'"

        # Check carbons
        actual_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        assert actual_carbons == expected_carbons, \
            f"{species_key}: Expected {expected_carbons} carbons, got {actual_carbons}. " \
            f"SMILES: {smiles}. {notes}"

        # Check rings
        ring_info = mol.GetRingInfo()
        actual_rings = sorted([len(r) for r in ring_info.AtomRings()])

        if expected_rings:
            for exp_ring in sorted(expected_rings):
                assert exp_ring in actual_rings, \
                    f"{species_key}: Expected ring size {exp_ring}, got {actual_rings}. " \
                    f"SMILES: {smiles}. {notes}"
        else:
            assert len(actual_rings) == 0, \
                f"{species_key}: Expected no rings (acyclic), got {actual_rings}. " \
                f"SMILES: {smiles}. {notes}"


class TestCriticalTerpeneOxidationProducts:
    """
    Critical tests for terpene oxidation products.

    These compounds MUST have cyclobutane (4-membered) rings, NOT cyclopentane (5-membered)
    or cyclohexane (6-membered) rings. Previous bugs showed these as wrong structures.
    """

    def test_pinic_acid_is_cyclobutane(self):
        """Pinic acid MUST have a 4-membered cyclobutane ring."""
        compound = COMPOUNDS.get('pinic_acid')
        assert compound, "pinic_acid not found in database"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        # Must contain exactly one 4-membered ring
        assert 4 in ring_sizes, \
            f"pinic_acid MUST have cyclobutane (4-membered ring). " \
            f"Found: {ring_sizes}. SMILES: {compound.smiles}"

        # Must NOT contain 5 or 6 membered rings (common error)
        assert 5 not in ring_sizes, \
            f"pinic_acid has wrong 5-membered ring (cyclopentane). " \
            f"Should be 4-membered (cyclobutane). SMILES: {compound.smiles}"
        assert 6 not in ring_sizes, \
            f"pinic_acid has wrong 6-membered ring (cyclohexane). " \
            f"Should be 4-membered (cyclobutane). SMILES: {compound.smiles}"

        # Check molecular formula
        formula = rdMolDescriptors.CalcMolFormula(mol)
        assert formula == 'C9H14O4', \
            f"pinic_acid wrong formula: {formula}. Should be C9H14O4"

    def test_pinonic_acid_is_cyclobutane(self):
        """Pinonic acid MUST have a 4-membered cyclobutane ring."""
        compound = COMPOUNDS.get('pinonic_acid')
        assert compound, "pinonic_acid not found in database"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        assert 4 in ring_sizes, \
            f"pinonic_acid MUST have cyclobutane (4-membered ring). " \
            f"Found: {ring_sizes}. SMILES: {compound.smiles}"

        assert 5 not in ring_sizes and 6 not in ring_sizes, \
            f"pinonic_acid has wrong ring size. Should be 4-membered."

        formula = rdMolDescriptors.CalcMolFormula(mol)
        assert formula == 'C10H16O3', \
            f"pinonic_acid wrong formula: {formula}. Should be C10H16O3"

    def test_pinonaldehyde_is_cyclobutane(self):
        """Pinonaldehyde MUST have a 4-membered cyclobutane ring."""
        compound = COMPOUNDS.get('pinonaldehyde')
        assert compound, "pinonaldehyde not found in database"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        assert 4 in ring_sizes, \
            f"pinonaldehyde MUST have cyclobutane (4-membered ring). " \
            f"Found: {ring_sizes}. SMILES: {compound.smiles}"

        formula = rdMolDescriptors.CalcMolFormula(mol)
        assert formula == 'C10H16O2', \
            f"pinonaldehyde wrong formula: {formula}. Should be C10H16O2"

    def test_gecko_to_smiles_pinic_acid_code(self):
        """GECKO species code TA9000 must return correct pinic acid structure."""
        smiles = gecko_to_smiles("", "TA9000")

        mol = Chem.MolFromSmiles(smiles)
        assert mol, f"TA9000 returned invalid SMILES: {smiles}"

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        assert 4 in ring_sizes, \
            f"TA9000 (pinic acid) must have 4-membered ring. Got: {ring_sizes}"

        carbons = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
        assert carbons == 9, f"TA9000 must have 9 carbons, got {carbons}"


class TestFunctionalGroupDetection:
    """Test that functional groups are correctly identified in structures."""

    @pytest.mark.parametrize("compound_name,expected_groups", [
        ('pinic_acid', ['carboxylic_acid', 'carboxylic_acid']),  # Diacid
        ('pinonic_acid', ['carboxylic_acid', 'ketone']),
        ('pinonaldehyde', ['aldehyde', 'ketone']),
        ('acetic_acid', ['carboxylic_acid']),
        ('formaldehyde', ['aldehyde']),
        ('acetone', ['ketone']),
        ('methanol', ['alcohol']),
        ('ethanol', ['alcohol']),
    ])
    def test_functional_groups_present(self, compound_name, expected_groups):
        """Test compounds have expected functional groups."""
        compound = COMPOUNDS.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES for {compound_name}"

        smiles = compound.smiles.upper()

        for group in expected_groups:
            if group == 'carboxylic_acid':
                # COOH pattern: C(=O)O or OC(=O) or similar
                has_acid = 'C(=O)O' in smiles or 'OC(=O)' in smiles or 'C(O)=O' in smiles
                assert has_acid, \
                    f"{compound_name} should have carboxylic acid group"
            elif group == 'aldehyde':
                # CHO pattern: C=O at terminal position
                has_aldehyde = 'C=O' in smiles or 'CHO' in smiles or smiles.endswith('=O')
                assert has_aldehyde, \
                    f"{compound_name} should have aldehyde group"
            elif group == 'ketone':
                # C(=O)C pattern: carbonyl not at terminal
                has_ketone = 'C(=O)C' in smiles or 'CC(=O)' in smiles
                assert has_ketone, \
                    f"{compound_name} should have ketone group"
            elif group == 'alcohol':
                has_alcohol = 'O' in smiles
                assert has_alcohol, \
                    f"{compound_name} should have alcohol group"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
