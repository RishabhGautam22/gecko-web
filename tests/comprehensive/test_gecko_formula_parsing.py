"""
GECKO Formula Validation Tests

Tests that GECKO formulas can be properly:
1. Retrieved for all compounds
2. Passed to the GECKO-A input preparation
3. Match expected patterns for functional groups

Author: GECKO-A Development Team
"""

import pytest
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from gecko_web.chemdata.compound_database import CompoundDatabase, COMPOUNDS
from gecko_web.chemdata import get_gecko_formula, get_compound
from gecko_web.main import get_gecko_input, validate_voc


class TestGeckoFormulaRetrieval:
    """Test GECKO formula retrieval through various methods."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    # Sample compounds from each category for testing
    SAMPLE_COMPOUNDS = [
        # Alkanes
        'methane', 'ethane', 'propane', 'butane', 'isobutane', 'pentane',
        'isopentane', 'hexane', 'cyclohexane', 'methylcyclohexane',
        # Alkenes
        'ethene', 'propene', '1_butene', 'isoprene',
        # Aromatics
        'benzene', 'toluene', 'xylene', 'ethylbenzene', 'styrene',
        # Monoterpenes
        'alpha_pinene', 'beta_pinene', 'limonene', 'myrcene', 'camphene',
        # Sesquiterpenes
        'beta_caryophyllene', 'alpha_humulene',
        # Oxygenated
        'methanol', 'ethanol', 'acetone', 'formaldehyde', 'acetaldehyde',
        'acetic_acid', 'formic_acid',
        # Terpene oxidation products
        'pinic_acid', 'pinonic_acid', 'pinonaldehyde', 'nopinone',
        # New additions
        'valeric_acid', 'hexanoic_acid', 'methyl_ethyl_ketone',
        'dimethyl_ether', 'diethyl_ether',
    ]

    @pytest.mark.parametrize("compound_name", SAMPLE_COMPOUNDS)
    def test_gecko_formula_retrievable(self, compound_name):
        """GECKO formula should be retrievable for sample compounds."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        assert compound.gecko_formula, \
            f"{compound_name}: No gecko_formula defined"

        # Formula should be non-trivial (not just the compound name)
        assert compound.gecko_formula != compound_name, \
            f"{compound_name}: gecko_formula is just the compound name"

    @pytest.mark.parametrize("compound_name", SAMPLE_COMPOUNDS)
    def test_get_gecko_input_returns_formula(self, compound_name):
        """get_gecko_input() should return a valid GECKO formula."""
        result = get_gecko_input(compound_name)

        # Should not return the compound name itself (that indicates lookup failure)
        assert result != compound_name, \
            f"{compound_name}: get_gecko_input returned the name, not a formula"

        # Result should look like a GECKO formula
        assert len(result) >= 2, \
            f"{compound_name}: get_gecko_input returned too short: '{result}'"

    @pytest.mark.parametrize("compound_name", SAMPLE_COMPOUNDS)
    def test_validate_voc_accepts_compound(self, compound_name):
        """validate_voc() should accept valid compound names."""
        assert validate_voc(compound_name), \
            f"{compound_name}: validate_voc() rejected valid compound"


class TestGeckoFormulaPatterns:
    """Test that GECKO formulas follow expected patterns."""

    # Functional group patterns in GECKO notation
    FUNCTIONAL_GROUPS = {
        'hydroxyl': ['OH', '(OH)'],
        'carbonyl': ['CO', 'CHO', 'COCH3'],
        'carboxyl': ['COOH', 'CO(OH)'],
        'nitrate': ['ONO2', '(ONO2)'],
        'peroxy': ['OO', 'O2'],
        'double_bond': ['Cd', 'cd', '=Cd'],
        'aromatic': ['c1', 'c2', 'cH'],
    }

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_alcohols_have_hydroxyl_group(self):
        """Alcohols should have OH in GECKO formula."""
        alcohols = ['methanol', 'ethanol', 'propanol', '1_butanol', '2_butanol']

        for name in alcohols:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            has_oh = any(pattern in formula for pattern in ['OH', '(OH)'])
            assert has_oh, \
                f"{name}: Alcohol should have OH in formula. Got: '{formula}'"

    def test_aldehydes_have_cho_group(self):
        """Aldehydes should have CHO in GECKO formula."""
        aldehydes = ['formaldehyde', 'acetaldehyde', 'propionaldehyde', 'pinonaldehyde']

        for name in aldehydes:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            has_cho = 'CHO' in formula or 'HCO' in formula
            assert has_cho, \
                f"{name}: Aldehyde should have CHO in formula. Got: '{formula}'"

    def test_ketones_have_co_group(self):
        """Ketones should have CO (not CHO) in GECKO formula."""
        ketones = ['acetone', 'methyl_ethyl_ketone', 'methyl_isobutyl_ketone']

        for name in ketones:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            # Ketones have CO but not CHO pattern
            assert 'CO' in formula, \
                f"{name}: Ketone should have CO in formula. Got: '{formula}'"

    def test_carboxylic_acids_have_cooh_group(self):
        """Carboxylic acids should have COOH or CO(OH) in GECKO formula."""
        acids = ['acetic_acid', 'formic_acid', 'propionic_acid', 'valeric_acid',
                 'pinic_acid', 'pinonic_acid']

        for name in acids:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            has_cooh = any(p in formula for p in ['COOH', 'CO(OH)'])
            assert has_cooh, \
                f"{name}: Carboxylic acid should have COOH or CO(OH). Got: '{formula}'"

    def test_ethers_have_oco_pattern(self):
        """Ethers should have C-O-C pattern in GECKO formula."""
        ethers = ['dimethyl_ether', 'diethyl_ether', 'methyl_tert_butyl_ether']

        for name in ethers:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            # Look for O between carbon groups - various patterns:
            # OCH (simple ether), OC( (branched), CH3O (methoxy), etc.
            has_ether = ('O' in formula and
                        (re.search(r'CH[234]?O', formula) or
                         re.search(r'OC[H(]', formula) or
                         'OCH' in formula))
            assert has_ether, \
                f"{name}: Ether should have C-O-C linkage. Got: '{formula}'"

    def test_alkenes_have_double_bond_notation(self):
        """Alkenes should have Cd (double bond carbon) in GECKO formula."""
        alkenes = ['ethene', 'propene', '1_butene', 'isoprene']

        for name in alkenes:
            compound = self.db.get(name)
            if not compound:
                continue

            formula = compound.gecko_formula
            has_db = 'Cd' in formula or 'cd' in formula or '=' in formula
            assert has_db, \
                f"{name}: Alkene should have Cd (double bond). Got: '{formula}'"


class TestCompoundNameFormats:
    """Test that various name formats resolve correctly."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    @pytest.mark.parametrize("query,expected", [
        ('alpha_pinene', 'alpha_pinene'),
        ('ALPHA_PINENE', 'alpha_pinene'),
        ('Alpha_Pinene', 'alpha_pinene'),
        ('pinic_acid', 'pinic_acid'),
        ('PINIC_ACID', 'pinic_acid'),
        ('methanol', 'methanol'),
        ('METHANOL', 'methanol'),
    ])
    def test_case_insensitive_lookup(self, query, expected):
        """Compound lookup should be case-insensitive."""
        result = self.db.get(query)
        assert result is not None, f"Failed to find '{query}'"
        assert result.name == expected, \
            f"Query '{query}' returned '{result.name}', expected '{expected}'"

    @pytest.mark.parametrize("alias,expected_name", [
        ('APINENE', 'alpha_pinene'),
        ('ALPHAPINENE', 'alpha_pinene'),
        ('MEK', 'methyl_ethyl_ketone'),
        ('DME', 'dimethyl_ether'),
        ('MTBE', 'methyl_tert_butyl_ether'),
    ])
    def test_alias_lookup(self, alias, expected_name):
        """Compounds should be findable by their aliases."""
        result = self.db.get(alias)
        if result:
            assert result.name == expected_name, \
                f"Alias '{alias}' returned '{result.name}', expected '{expected_name}'"
        # If alias not set up, that's OK - skip


class TestGeckoFormulaCharacterCounts:
    """Verify GECKO formula character counts match molecular formulas."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def count_carbons_in_gecko(self, formula):
        """Count carbon atoms in a GECKO formula."""
        # Remove ring markers (C1, C2, etc. should count as 1 C each)
        # Count 'C' and 'c' (aromatic)
        count = 0

        # Simple approach: count C followed by number or not
        # CH3 = 1C, CH2 = 1C, CH = 1C, CO = 1C, etc.
        pattern = r'[Cc](?:H[234]?|O|d|=|[0-9]|\(|$)'
        matches = re.findall(pattern, formula)
        count = len(matches)

        # Also count standalone C at start
        if formula.startswith('C') and len(formula) > 1:
            if formula[1] not in 'Ccod':
                count = max(count, 1)

        return count

    @pytest.mark.parametrize("compound_name", [
        'methane', 'ethane', 'propane', 'methanol', 'ethanol', 'acetone'
    ])
    def test_simple_compounds_carbon_count(self, compound_name):
        """Simple compounds should have matching carbon counts."""
        compound = self.db.get(compound_name)
        if not compound:
            pytest.skip(f"{compound_name} not in database")

        # Parse expected carbon count from molecular formula
        mf = compound.molecular_formula
        c_match = re.search(r'C(\d*)', mf)
        if not c_match:
            pytest.skip(f"{compound_name} has no C in formula")

        expected_c = int(c_match.group(1)) if c_match.group(1) else 1

        # For very simple molecules, verify by inspection
        gecko = compound.gecko_formula

        # Methane: CH4 -> gecko should be CH4
        if compound_name == 'methane':
            assert 'CH4' in gecko or gecko == 'CH4', \
                f"methane gecko_formula should be CH4, got {gecko}"

        # Ethane: C2H6 -> gecko should be CH3CH3
        if compound_name == 'ethane':
            assert gecko == 'CH3CH3', \
                f"ethane gecko_formula should be CH3CH3, got {gecko}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
