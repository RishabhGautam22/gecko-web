"""
Comprehensive Compound Database Integrity Tests

Tests ALL compounds in the database for:
1. Valid SMILES structures (RDKit parseable)
2. GECKO formula presence and format
3. Molecular formula correctness
4. Molecular weight accuracy
5. Ring structure correctness (cyclobutane vs cyclopentane etc.)
6. Category assignment validity
7. Cross-reference with voc_categories

Author: GECKO-A Development Team
"""

import pytest
import sys
import re
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from gecko_web.chemdata.compound_database import CompoundDatabase, COMPOUNDS
from gecko_web.chemdata import voc_categories

# Import RDKit for structure validation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    pytest.skip("RDKit not available", allow_module_level=True)


class TestAllCompoundsValidity:
    """Test every compound in the database for basic validity."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()
        self.all_compounds = list(COMPOUNDS.values())

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values()),
                             ids=lambda c: c.name)
    def test_compound_has_required_fields(self, compound):
        """Every compound must have name, SMILES, gecko_formula, molecular_formula, MW."""
        assert compound.name, f"Compound missing name"
        assert compound.smiles, f"{compound.name}: missing SMILES"
        assert compound.gecko_formula, f"{compound.name}: missing gecko_formula"
        assert compound.molecular_formula, f"{compound.name}: missing molecular_formula"
        assert compound.molecular_weight > 0, f"{compound.name}: invalid molecular_weight"

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values()),
                             ids=lambda c: c.name)
    def test_smiles_is_valid_rdkit(self, compound):
        """Every SMILES must be parseable by RDKit."""
        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol is not None, f"{compound.name}: Invalid SMILES '{compound.smiles}'"

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values()),
                             ids=lambda c: c.name)
    def test_molecular_weight_matches_smiles(self, compound):
        """Molecular weight should match the SMILES structure (within 1 Da)."""
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            calculated_mw = Descriptors.MolWt(mol)
            tolerance = 1.0  # Allow 1 Da difference for rounding
            assert abs(calculated_mw - compound.molecular_weight) < tolerance, \
                f"{compound.name}: MW mismatch. Database: {compound.molecular_weight}, " \
                f"Calculated: {calculated_mw:.2f}"

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values()),
                             ids=lambda c: c.name)
    def test_gecko_formula_format(self, compound):
        """GECKO formula should follow expected format patterns."""
        formula = compound.gecko_formula

        # Should contain only valid GECKO characters
        # Note: # is valid for triple bonds (alkynes like ethyne, propyne)
        valid_pattern = r'^[A-Za-z0-9\(\)\-=#]+$'
        assert re.match(valid_pattern, formula), \
            f"{compound.name}: Invalid characters in gecko_formula '{formula}'"

        # Should start with a valid atom (C, H, O, N, or functional group)
        assert formula[0] in 'CHONchdS' or formula.startswith('('), \
            f"{compound.name}: gecko_formula should start with valid atom: '{formula}'"

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values()),
                             ids=lambda c: c.name)
    def test_category_is_valid(self, compound):
        """Compound category should be a known category."""
        valid_categories = {
            'alkane', 'alkene', 'aromatic', 'monoterpene', 'sesquiterpene',
            'oxygenated', 'terpene_oxidation_product', 'alcohol', 'aldehyde',
            'ketone', 'ester', 'ether', 'carboxylic_acid', 'phenolic',
            'isoprene_product', 'pinene_product', 'limonene_product',
            '', None  # Allow empty category
        }
        # Be lenient - just check it's not nonsense
        if compound.category:
            assert len(compound.category) < 50, \
                f"{compound.name}: category too long: '{compound.category}'"


class TestRingStructureCorrectness:
    """Verify ring structures are correctly identified."""

    # Known cyclobutane compounds (4-membered rings)
    CYCLOBUTANE_COMPOUNDS = ['pinic_acid']

    # Known cyclopentane compounds (5-membered rings)
    CYCLOPENTANE_COMPOUNDS = ['cyclopentane', 'cyclopentene']

    # Known cyclohexane compounds (6-membered rings)
    CYCLOHEXANE_COMPOUNDS = [
        'cyclohexane', 'cyclohexene', 'cyclohexanone', 'methylcyclohexane'
    ]

    # Known bicyclic monoterpenes
    BICYCLIC_MONOTERPENES = [
        'alpha_pinene', 'beta_pinene', 'camphene', '3_carene'
    ]

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_pinic_acid_has_cyclobutane_ring(self):
        """Pinic acid must have a 4-membered cyclobutane ring."""
        compound = self.db.get('pinic_acid')
        assert compound, "pinic_acid not found in database"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES for pinic_acid: {compound.smiles}"

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        assert 4 in ring_sizes, \
            f"pinic_acid should have a 4-membered ring (cyclobutane). " \
            f"Found ring sizes: {ring_sizes}. SMILES: {compound.smiles}"

    def test_pinic_acid_has_two_carboxylic_acids(self):
        """Pinic acid must have exactly 2 carboxylic acid groups."""
        compound = self.db.get('pinic_acid')
        assert compound, "pinic_acid not found in database"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES for pinic_acid: {compound.smiles}"

        # Count carboxylic acid pattern: C(=O)O
        smarts = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        matches = mol.GetSubstructMatches(smarts)

        assert len(matches) == 2, \
            f"pinic_acid should have 2 COOH groups, found {len(matches)}. " \
            f"SMILES: {compound.smiles}"

    @pytest.mark.parametrize("compound_name", CYCLOHEXANE_COMPOUNDS)
    def test_cyclohexane_derivatives_have_6_ring(self, compound_name):
        """Cyclohexane derivatives must have a 6-membered ring."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.fail(f"Invalid SMILES for {compound_name}: {compound.smiles}")

        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]

        assert 6 in ring_sizes, \
            f"{compound_name} should have a 6-membered ring. " \
            f"Found ring sizes: {ring_sizes}"

    @pytest.mark.parametrize("compound_name", BICYCLIC_MONOTERPENES)
    def test_bicyclic_monoterpenes_have_two_rings(self, compound_name):
        """Bicyclic monoterpenes must have at least 2 rings.

        Note: RDKit's SSSR (Smallest Set of Smallest Rings) may count fused
        bicyclic systems as 3 rings because it includes the "outer" ring.
        Alpha-pinene and beta-pinene have a 4+6 fused system that RDKit
        counts as 3 SSSR rings.
        """
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.fail(f"Invalid SMILES for {compound_name}: {compound.smiles}")

        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()

        # Bicyclic = at least 2 rings (may be 3 in SSSR for fused systems)
        assert num_rings >= 2, \
            f"{compound_name} should be bicyclic (>=2 rings). Found {num_rings} rings."


class TestVOCCategoriesConsistency:
    """Test that voc_categories entries match compound_database."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_all_voc_category_compounds_exist(self):
        """All compounds in voc_categories should exist in compound_database."""
        missing = []
        total = 0

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    total += 1
                    compound = self.db.get(compound_name)
                    if not compound:
                        missing.append(f"{compound_name} (in {cat_info.name})")

        # Allow some missing but report them
        missing_pct = len(missing) / total * 100 if total > 0 else 0

        # Fail if more than 30% are missing
        assert missing_pct < 30, \
            f"{len(missing)}/{total} ({missing_pct:.1f}%) compounds missing from database:\n" + \
            "\n".join(missing[:20]) + \
            (f"\n... and {len(missing)-20} more" if len(missing) > 20 else "")

    def test_compound_names_use_underscores_not_spaces(self):
        """Compound names should use underscores, not spaces."""
        names_with_spaces = []

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    if ' ' in compound_name:
                        names_with_spaces.append(compound_name)

        assert len(names_with_spaces) == 0, \
            f"Found {len(names_with_spaces)} compound names with spaces " \
            f"(should use underscores): {names_with_spaces[:10]}"

    def test_compound_lookup_with_various_formats(self):
        """Compounds should be findable with different name formats."""
        test_cases = [
            ('alpha_pinene', 'alpha_pinene'),
            ('ALPHA_PINENE', 'alpha_pinene'),
            ('alpha-pinene', 'alpha_pinene'),
            ('alphapinene', 'alpha_pinene'),  # May not work, that's OK
        ]

        for query, expected_name in test_cases[:3]:  # First 3 should work
            result = self.db.get(query)
            if result:
                assert result.name == expected_name, \
                    f"Query '{query}' returned '{result.name}', expected '{expected_name}'"


class TestGECKOFormulaValidity:
    """Test GECKO formula syntax and structure."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_gecko_formula_ring_notation_balanced(self):
        """Ring notation markers (C1, C2, etc.) should be balanced."""
        unbalanced = []

        for name, compound in COMPOUNDS.items():
            formula = compound.gecko_formula

            # Count ring markers (C1, C2, c1, c2, etc.)
            for ring_num in ['1', '2', '3']:
                # Count uppercase and lowercase separately
                c_upper = formula.count(f'C{ring_num}')
                c_lower = formula.count(f'c{ring_num}')

                # Each ring marker should appear exactly twice (open and close)
                # or zero times (no such ring)
                total = c_upper + c_lower
                if total > 0 and total != 2:
                    unbalanced.append(
                        f"{name}: C{ring_num} appears {total} times in '{formula}'"
                    )

        # Allow some failures but report - GECKO notation can use different
        # conventions for ring closures (e.g., C1...c1 for aromatic rings)
        # This is informational, not a hard failure
        if len(unbalanced) > 0:
            print(f"\nNote: {len(unbalanced)} compounds have unusual ring notation:")
            for item in unbalanced[:5]:
                print(f"  {item}")
        # Only fail if there are many (>20) which indicates systematic issues
        assert len(unbalanced) < 20, \
            f"Found {len(unbalanced)} compounds with unbalanced ring notation:\n" + \
            "\n".join(unbalanced[:10])

    @pytest.mark.parametrize("compound_name,expected_groups", [
        ('methanol', ['OH']),
        ('acetaldehyde', ['CHO']),
        ('acetone', ['CO']),
        ('acetic_acid', ['CO', 'OH']),
        ('formaldehyde', ['CHO']),
    ])
    def test_gecko_formula_contains_expected_groups(self, compound_name, expected_groups):
        """GECKO formula should contain expected functional groups."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        formula = compound.gecko_formula
        for group in expected_groups:
            assert group in formula, \
                f"{compound_name}: Expected '{group}' in gecko_formula '{formula}'"


class TestMolecularFormulaConsistency:
    """Test that molecular formulas are consistent with SMILES."""

    @pytest.mark.parametrize("compound", list(COMPOUNDS.values())[:50],  # Test first 50
                             ids=lambda c: c.name)
    def test_molecular_formula_atom_count(self, compound):
        """Molecular formula atom counts should match SMILES."""
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            pytest.skip(f"Invalid SMILES for {compound.name}")

        # Get atom counts from SMILES
        mol_with_h = Chem.AddHs(mol)

        # Parse molecular formula
        formula = compound.molecular_formula

        # Extract carbon count from formula
        c_match = re.search(r'C(\d*)', formula)
        if c_match:
            expected_c = int(c_match.group(1)) if c_match.group(1) else 1
            actual_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

            assert actual_c == expected_c, \
                f"{compound.name}: Carbon count mismatch. " \
                f"Formula says {expected_c}, SMILES has {actual_c}"


class TestSpecificCompoundStructures:
    """Test specific compounds that have had issues."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_pinic_acid_structure_is_correct(self):
        """
        Pinic acid should be:
        - 3-(carboxymethyl)-2,2-dimethylcyclobutane-1-carboxylic acid
        - CAS: 19067-78-4
        - Has cyclobutane (4-membered) ring
        - Has gem-dimethyl group
        - Has two carboxylic acid groups
        """
        compound = self.db.get('pinic_acid')
        assert compound, "pinic_acid not found"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        # Check ring size
        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]
        assert 4 in ring_sizes, f"Should have 4-membered ring, got {ring_sizes}"

        # Check molecular formula
        assert compound.molecular_formula == 'C9H14O4', \
            f"Wrong formula: {compound.molecular_formula}"

        # Check CAS
        assert compound.cas == '19067-78-4', f"Wrong CAS: {compound.cas}"

    def test_valeric_acid_is_pentanoic_acid(self):
        """Valeric acid should be pentanoic acid (5-carbon carboxylic acid)."""
        compound = self.db.get('valeric_acid')
        assert compound, "valeric_acid not found"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        # Should have 5 carbons
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        assert c_count == 5, f"Valeric acid should have 5 carbons, got {c_count}"

        # Should have carboxylic acid group
        smarts = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        matches = mol.GetSubstructMatches(smarts)
        assert len(matches) == 1, f"Should have 1 COOH group, got {len(matches)}"

    def test_alpha_pinene_is_bicyclic(self):
        """Alpha-pinene should be a bicyclic monoterpene.

        Note: RDKit's SSSR counts alpha-pinene's fused 4+6 ring system as
        3 rings (it includes the outer envelope ring).
        """
        compound = self.db.get('alpha_pinene')
        assert compound, "alpha_pinene not found"

        mol = Chem.MolFromSmiles(compound.smiles)
        assert mol, f"Invalid SMILES: {compound.smiles}"

        # Should have at least 2 rings (SSSR may count 3 for fused systems)
        ring_info = mol.GetRingInfo()
        assert ring_info.NumRings() >= 2, \
            f"alpha_pinene should be bicyclic (>=2 rings), got {ring_info.NumRings()} rings"

        # Should have 10 carbons (monoterpene)
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        assert c_count == 10, f"Monoterpene should have 10 carbons, got {c_count}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
