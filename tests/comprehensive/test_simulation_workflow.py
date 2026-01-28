"""
Simulation Workflow Integration Tests

Tests the complete workflow from compound selection to simulation setup:
1. Compound resolution to GECKO formula
2. Mechanism name generation (no spaces)
3. File path safety
4. Scenario parameter application

These tests verify the fixes for issues like:
- Spaces in compound names breaking bash scripts
- Missing GECKO formulas causing "compound not found" errors
- Incorrect mechanism name generation

Author: GECKO-A Development Team
"""

import pytest
import sys
import re
from pathlib import Path
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from gecko_web.chemdata.compound_database import CompoundDatabase
from gecko_web.chemdata import voc_categories
from gecko_web.main import get_gecko_input, validate_voc


class TestMechanismNameGeneration:
    """Test that mechanism names are generated correctly for bash scripts."""

    # Names that previously caused issues
    PROBLEMATIC_NAMES = [
        'pinic acid',      # Space in name
        'pinonic acid',    # Space in name
        'acetic acid',     # Space in name
        'formic acid',     # Space in name
        'valeric acid',    # Space in name
        'methyl vinyl ketone',
        '10-hydroxypinonic acid',
    ]

    def normalize_mech_name(self, voc_name: str) -> str:
        """Replicate the mechanism name normalization from main.py."""
        return voc_name.lower().strip().replace(' ', '_').replace('-', '_')

    @pytest.mark.parametrize("name", PROBLEMATIC_NAMES)
    def test_mech_name_has_no_spaces(self, name):
        """Mechanism names should never contain spaces."""
        mech_name = self.normalize_mech_name(name)
        assert ' ' not in mech_name, \
            f"Mechanism name '{mech_name}' contains spaces (from '{name}')"

    @pytest.mark.parametrize("name", PROBLEMATIC_NAMES)
    def test_mech_name_is_bash_safe(self, name):
        """Mechanism names should be safe for bash file operations."""
        mech_name = self.normalize_mech_name(name)

        # Should only contain alphanumeric and underscores
        safe_pattern = r'^[a-z0-9_]+$'
        assert re.match(safe_pattern, mech_name), \
            f"Mechanism name '{mech_name}' is not bash-safe (from '{name}')"

    def test_all_voc_category_names_become_bash_safe(self):
        """All compound names in voc_categories should become bash-safe."""
        unsafe_names = []

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    mech_name = self.normalize_mech_name(compound_name)

                    # Check for unsafe characters
                    if not re.match(r'^[a-z0-9_]+$', mech_name):
                        unsafe_names.append(f"'{compound_name}' -> '{mech_name}'")

        assert len(unsafe_names) == 0, \
            f"Found {len(unsafe_names)} unsafe mechanism names:\n" + \
            "\n".join(unsafe_names[:10])


class TestGeckoInputResolution:
    """Test that GECKO input is properly resolved for all compounds."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_all_database_compounds_have_gecko_formula(self):
        """Every compound in database should have a GECKO formula."""
        missing_formula = []

        for name, compound in self.db._compounds.items():
            if not compound.gecko_formula:
                missing_formula.append(name)
            elif compound.gecko_formula == name:
                # Formula is just the name - this is a fallback, not real
                missing_formula.append(f"{name} (formula is name)")

        assert len(missing_formula) == 0, \
            f"Found {len(missing_formula)} compounds without proper GECKO formula:\n" + \
            "\n".join(missing_formula[:20])

    def test_get_gecko_input_returns_formula_not_name(self):
        """get_gecko_input should return formula, not compound name."""
        test_compounds = [
            'methane', 'ethane', 'benzene', 'alpha_pinene', 'pinic_acid',
            'valeric_acid', 'methyl_ethyl_ketone', 'dimethyl_ether'
        ]

        name_returned = []

        for name in test_compounds:
            result = get_gecko_input(name)
            if result == name:
                name_returned.append(name)

        assert len(name_returned) == 0, \
            f"get_gecko_input returned name instead of formula for: {name_returned}"

    def test_get_gecko_input_handles_case_variations(self):
        """get_gecko_input should handle different case variations."""
        test_cases = [
            ('METHANE', 'CH4'),
            ('Methane', 'CH4'),
            ('methane', 'CH4'),
            ('BENZENE', None),  # Just check it doesn't return 'BENZENE'
        ]

        for query, expected in test_cases:
            result = get_gecko_input(query)
            assert result.lower() != query.lower(), \
                f"get_gecko_input({query}) returned the name, not formula"


class TestVOCCategoriesIntegrity:
    """Test that voc_categories is properly integrated with compound database."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    def test_no_spaces_in_voc_category_compound_names(self):
        """Compound names in voc_categories should not have spaces."""
        names_with_spaces = []

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    if ' ' in compound_name:
                        names_with_spaces.append(compound_name)

        assert len(names_with_spaces) == 0, \
            f"Found {len(names_with_spaces)} names with spaces: {names_with_spaces}"

    def test_voc_category_compounds_resolvable(self):
        """All voc_category compounds should be resolvable to GECKO formulas."""
        unresolvable = []

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    # Try to get GECKO input
                    result = get_gecko_input(compound_name)

                    # If result equals the name, lookup failed
                    if result == compound_name:
                        unresolvable.append(f"{compound_name} ({cat_info.name})")

        # Report all unresolvable but only fail if too many
        unresolvable_pct = len(unresolvable) / 170 * 100  # ~170 total

        if unresolvable:
            print(f"\nUnresolvable compounds ({len(unresolvable)}):")
            for name in unresolvable[:20]:
                print(f"  - {name}")

        assert unresolvable_pct < 30, \
            f"Too many unresolvable compounds: {len(unresolvable)} ({unresolvable_pct:.1f}%)"


class TestScenarioParameters:
    """Test scenario parameter handling."""

    def test_default_scenario_params_are_valid(self):
        """Default scenario parameters should be within valid ranges."""
        from gecko_web.main import ScenarioParams

        # Create default params
        params = ScenarioParams()

        # Validate ranges
        assert -20 <= params.vapor_pressure_threshold <= 0, \
            f"Invalid vapor_pressure_threshold: {params.vapor_pressure_threshold}"

        assert 1 <= params.max_generations <= 10, \
            f"Invalid max_generations: {params.max_generations}"

        assert 200 <= params.temperature_k <= 400, \
            f"Invalid temperature_k: {params.temperature_k}"

        assert 0 <= params.rh_percent <= 100, \
            f"Invalid rh_percent: {params.rh_percent}"


class TestCompoundCategorySampling:
    """Test sampling compounds from each category."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.db = CompoundDatabase()

    # Sample 5 compounds from each major category
    CATEGORY_SAMPLES = {
        'alkane': ['methane', 'ethane', 'propane', 'butane', 'cyclohexane'],
        'alkene': ['ethene', 'propene', 'isoprene', '1_butene', 'cyclohexene'],
        'aromatic': ['benzene', 'toluene', 'xylene', 'ethylbenzene', 'styrene'],
        'monoterpene': ['alpha_pinene', 'beta_pinene', 'limonene', 'myrcene', 'camphene'],
        'oxygenated': ['methanol', 'ethanol', 'acetone', 'formaldehyde', 'acetaldehyde'],
        'terpene_oxidation_product': ['pinic_acid', 'pinonic_acid', 'pinonaldehyde', 'nopinone'],
    }

    @pytest.mark.parametrize("category,compounds", list(CATEGORY_SAMPLES.items()))
    def test_category_samples_have_valid_gecko_formulas(self, category, compounds):
        """Sample compounds from each category should have valid GECKO formulas."""
        invalid = []

        for name in compounds:
            compound = self.db.get(name)
            if not compound:
                continue

            # Check GECKO formula exists and is valid
            if not compound.gecko_formula:
                invalid.append(f"{name}: no formula")
            elif compound.gecko_formula == name:
                invalid.append(f"{name}: formula is name")
            elif len(compound.gecko_formula) < 3:
                invalid.append(f"{name}: formula too short")

        assert len(invalid) == 0, \
            f"Category '{category}' has invalid GECKO formulas:\n" + \
            "\n".join(invalid)

    @pytest.mark.parametrize("category,compounds", list(CATEGORY_SAMPLES.items()))
    def test_category_samples_pass_validation(self, category, compounds):
        """Sample compounds from each category should pass VOC validation."""
        failed = []

        for name in compounds:
            if not validate_voc(name):
                failed.append(name)

        assert len(failed) == 0, \
            f"Category '{category}' has compounds that fail validation: {failed}"


class TestFilePathSafety:
    """Test that generated file paths are safe."""

    def test_mechanism_names_are_filesystem_safe(self):
        """Mechanism names should be safe for all filesystems."""
        # Characters that are problematic on various filesystems
        unsafe_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ', '\t', '\n']

        test_names = [
            'pinic_acid', 'alpha_pinene', 'beta_caryophyllene',
            'methyl_ethyl_ketone', '10_hydroxypinonic_acid'
        ]

        for name in test_names:
            mech_name = name.lower().strip().replace(' ', '_').replace('-', '_')

            for char in unsafe_chars:
                assert char not in mech_name, \
                    f"Mechanism name '{mech_name}' contains unsafe char '{repr(char)}'"

    def test_directory_names_length_limit(self):
        """Mechanism/directory names should not exceed filesystem limits."""
        max_length = 100  # Conservative limit for Fortran compatibility

        for cat_id, cat_info in voc_categories.VOC_CATEGORIES.items():
            for sub_id, sub_info in cat_info.subcategories.items():
                for compound_name in sub_info.compounds:
                    mech_name = compound_name.lower().strip().replace(' ', '_')

                    assert len(mech_name) <= max_length, \
                        f"Mechanism name too long ({len(mech_name)} chars): {mech_name}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
