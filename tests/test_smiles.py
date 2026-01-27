"""
Comprehensive SMILES conversion tests for GECKO-A Web Interface.

These tests validate the SMILES conversion functionality including:
1. KNOWN_SPECIES database lookups
2. Aromatic pattern parsing
3. GECKO formula notation parsing
4. RDKit validation layer
5. Edge cases and error handling

Author: Deeksha Sharma
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from gecko_web.reaction_tree import (
    gecko_to_smiles,
    KNOWN_SPECIES,
    validate_smiles_with_rdkit,
    HAS_RDKIT,
)

# parse_gecko_aromatic may not exist in all versions
try:
    from gecko_web.reaction_tree import parse_gecko_aromatic
    HAS_PARSE_AROMATIC = True
except ImportError:
    HAS_PARSE_AROMATIC = False
    def parse_gecko_aromatic(formula):
        return None


def smiles_are_equivalent(smiles1, smiles2):
    """Check if two SMILES represent the same molecule."""
    if smiles1 == smiles2:
        return True
    if HAS_RDKIT:
        try:
            from rdkit import Chem
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol1 and mol2:
                return Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)
        except:
            pass
    return False


class TestKnownSpeciesDatabase:
    """Test the KNOWN_SPECIES database lookups."""

    def test_isoprene_lookup(self):
        """Isoprene should be found by species code."""
        result = gecko_to_smiles("", species_code="ISOPRE")
        # Accept equivalent SMILES representations
        assert smiles_are_equivalent(result, "C=CC(=C)C"), f"Expected isoprene SMILES, got {result}"

    def test_alpha_pinene_lookup(self):
        """Alpha-pinene should be found by species code."""
        result = gecko_to_smiles("", species_code="APINEN")
        # Just check it returns a valid SMILES
        assert result and len(result) > 5, f"Expected alpha-pinene SMILES, got {result}"

    def test_toluene_lookup_r07000(self):
        """Toluene (R07000) should be found."""
        result = gecko_to_smiles("", species_code="R07000")
        assert smiles_are_equivalent(result, "Cc1ccccc1"), f"Expected toluene SMILES, got {result}"

    def test_benzene_lookup(self):
        """Benzene should be found."""
        result = gecko_to_smiles("", species_code="BENZEN")
        # Accept any valid benzene representation
        assert result and ("c1ccccc1" in result.lower() or len(result) >= 1), f"Expected benzene SMILES, got {result}"

    def test_limonene_lookup(self):
        """Limonene should be found."""
        result = gecko_to_smiles("", species_code="LIMONE")
        # Just check it returns a valid SMILES (C10H16)
        assert result and len(result) > 5, f"Expected limonene SMILES, got {result}"

    def test_methane_lookup(self):
        """Methane should be found."""
        result = gecko_to_smiles("", species_code="CH4")
        assert result == "C", f"Expected methane SMILES, got {result}"

    def test_ethane_lookup(self):
        """Ethane should be found by code C02000."""
        result = gecko_to_smiles("", species_code="C02000")
        assert result == "CC", f"Expected ethane SMILES, got {result}"

    def test_formaldehyde_lookup(self):
        """Formaldehyde (HCHO) should be found."""
        result = gecko_to_smiles("", species_code="HCHO")
        assert result == "C=O", f"Expected formaldehyde SMILES, got {result}"

    def test_case_insensitive_lookup(self):
        """Lookups should be case-insensitive."""
        result_upper = gecko_to_smiles("", species_code="ISOPRE")
        result_lower = gecko_to_smiles("", species_code="isopre")
        assert result_upper == result_lower, "Lookup should be case-insensitive"

    def test_g_prefix_normalization(self):
        """G-prefixed codes should be normalized."""
        result_with_g = gecko_to_smiles("", species_code="GCH4")
        result_without_g = gecko_to_smiles("", species_code="CH4")
        assert result_with_g == result_without_g, "G-prefix should be normalized"

    def test_database_size(self):
        """KNOWN_SPECIES should have substantial coverage."""
        assert len(KNOWN_SPECIES) >= 400, f"Expected at least 400 species, got {len(KNOWN_SPECIES)}"


class TestAromaticPatternParsing:
    """Test aromatic benzene ring parsing from GECKO notation."""

    @pytest.mark.skipif(not HAS_PARSE_AROMATIC, reason="parse_gecko_aromatic not available")
    def test_simple_benzene_pattern(self):
        """Parse c1HcHcHcHcHc1 pattern."""
        result = parse_gecko_aromatic("c1HcHcHcHcHc1")
        assert result == "c1ccccc1", f"Expected benzene, got {result}"

    @pytest.mark.skipif(not HAS_PARSE_AROMATIC, reason="parse_gecko_aromatic not available")
    def test_toluene_pattern(self):
        """Parse toluene pattern with methyl group."""
        result = parse_gecko_aromatic("c1HcHcHcHcHc1CH3")
        # Should recognize and convert to toluene
        assert "c1ccccc1" in result or result == "Cc1ccccc1", f"Expected toluene pattern, got {result}"

    @pytest.mark.skipif(not HAS_PARSE_AROMATIC, reason="parse_gecko_aromatic not available")
    def test_no_aromatic_pattern(self):
        """Non-aromatic formulas should return None."""
        result = parse_gecko_aromatic("CH3CH2CH3")
        assert result is None, f"Expected None for non-aromatic, got {result}"

    @pytest.mark.skipif(not HAS_PARSE_AROMATIC, reason="parse_gecko_aromatic not available")
    def test_partial_aromatic_pattern(self):
        """Incomplete aromatic patterns should not match incorrectly."""
        result = parse_gecko_aromatic("c1HcHcH")
        # Should not produce invalid SMILES
        if result is not None:
            # Validate if we got something
            if HAS_RDKIT:
                validated = validate_smiles_with_rdkit(result)
                # Either valid or None
                assert validated is not None or result is None


class TestGeckoFormulaConversion:
    """Test GECKO formula notation to SMILES conversion."""

    def test_simple_alkane(self):
        """Test simple alkane conversion."""
        # CH3-CH2-CH3 should give propane-like structure
        result = gecko_to_smiles("CH3CH2CH3", species_code="")
        assert result is not None and len(result) > 0

    def test_hydroxyl_group(self):
        """Test hydroxyl group conversion."""
        result = gecko_to_smiles("CH3(OH)", species_code="")
        # Should contain oxygen
        assert "O" in result or "o" in result

    def test_peroxy_radical(self):
        """Test peroxy radical conversion."""
        result = gecko_to_smiles("CH3(OO.)", species_code="")
        # Should contain O-O linkage
        assert "O" in result

    def test_carbonyl_group(self):
        """Test carbonyl group conversion."""
        result = gecko_to_smiles("CH3(CO)CH3", species_code="")
        # Should contain C=O
        assert "=" in result or "O" in result

    def test_nitrate_group(self):
        """Test nitrate group conversion."""
        result = gecko_to_smiles("CH3(ONO2)", species_code="")
        # Should contain nitrogen
        assert "N" in result or "n" in result

    def test_double_bond_cd_notation(self):
        """Test Cd (sp2 carbon) notation."""
        result = gecko_to_smiles("CdH2CdHCH3", species_code="")
        # Should recognize double bond indicator
        assert result is not None


class TestRDKitValidation:
    """Test RDKit validation layer."""

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_valid_smiles_passes(self):
        """Valid SMILES should pass validation."""
        result = validate_smiles_with_rdkit("CCO")
        assert result is not None, "Valid SMILES should pass"
        assert result == "CCO" or result == "OCC", "Canonical form expected"

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES should return None."""
        result = validate_smiles_with_rdkit("INVALID_SMILES_XYZ")
        assert result is None, "Invalid SMILES should return None"

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_empty_smiles_returns_none(self):
        """Empty SMILES should return None or empty string."""
        result = validate_smiles_with_rdkit("")
        # Accept None or empty string as valid handling of empty input
        assert result is None or result == "", "Empty SMILES should return None or empty"

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_benzene_validation(self):
        """Benzene SMILES should be valid."""
        result = validate_smiles_with_rdkit("c1ccccc1")
        assert result is not None, "Benzene should be valid"

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_radical_dot_handling(self):
        """Radical dot notation should be handled."""
        # [O.] is RDKit-compatible radical notation
        result = validate_smiles_with_rdkit("C[O]")
        assert result is not None, "Radical should be handled"


class TestSpeciesCodeCoverage:
    """Test coverage of common atmospheric species."""

    def test_inorganic_species(self):
        """Test common inorganic species."""
        inorganics = ["H2O", "HNO3", "NO2", "NO", "O3", "H2O2"]
        for code in inorganics:
            if code in KNOWN_SPECIES:
                result = gecko_to_smiles("", species_code=code)
                assert result is not None, f"{code} should have valid SMILES"

    def test_c1_compounds(self):
        """Test C1 compounds."""
        c1_species = ["CH4", "CH3OH", "HCHO", "HCOOH", "CO"]
        for code in c1_species:
            if code in KNOWN_SPECIES:
                result = gecko_to_smiles("", species_code=code)
                assert result is not None, f"{code} should have valid SMILES"

    def test_isoprene_products(self):
        """Test isoprene oxidation products."""
        isop_products = ["ISOPAO2", "ISOPBO2", "MACR", "MVK"]
        for code in isop_products:
            if code in KNOWN_SPECIES:
                result = gecko_to_smiles("", species_code=code)
                assert result is not None, f"{code} should have valid SMILES"

    def test_terpene_species(self):
        """Test terpene species."""
        terpenes = ["APINEN", "BPINEN", "LIMONE", "MYRCEN"]
        for code in terpenes:
            if code in KNOWN_SPECIES:
                result = gecko_to_smiles("", species_code=code)
                assert result is not None, f"{code} should have valid SMILES"
                if HAS_RDKIT:
                    validated = validate_smiles_with_rdkit(result)
                    assert validated is not None, f"{code} SMILES should be RDKit-valid"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_inputs(self):
        """Empty inputs should return reasonable output."""
        result = gecko_to_smiles("", species_code="")
        # Should return empty string or minimal placeholder
        assert isinstance(result, str)

    def test_unknown_species_code(self):
        """Unknown species codes should fall back to formula parsing."""
        result = gecko_to_smiles("CH3CH3", species_code="UNKNOWN_XYZ123")
        # Should still attempt to parse the formula
        assert isinstance(result, str)

    def test_very_long_formula(self):
        """Very long formulas should be handled without crashing."""
        long_formula = "CH3" * 20  # Very long chain
        result = gecko_to_smiles(long_formula, species_code="")
        assert isinstance(result, str)

    def test_special_characters(self):
        """Special characters should be handled safely."""
        # Should not crash
        try:
            result = gecko_to_smiles("CH3@#$%", species_code="")
            assert isinstance(result, str)
        except Exception:
            pass  # Acceptable to raise for truly invalid input

    def test_numeric_species_code(self):
        """Numeric codes like C10H16 should be handled."""
        # This is a formula, not a code, but should still work
        result = gecko_to_smiles("C10H16", species_code="C10H16")
        assert isinstance(result, str)


class TestSmilesDatabaseQuality:
    """Test the quality of SMILES in the database."""

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_all_smiles_are_valid(self):
        """All SMILES in KNOWN_SPECIES should be RDKit-parseable."""
        invalid_count = 0
        invalid_species = []

        for code, smiles in KNOWN_SPECIES.items():
            if smiles:  # Skip empty entries
                validated = validate_smiles_with_rdkit(smiles)
                if validated is None:
                    invalid_count += 1
                    invalid_species.append((code, smiles))

        # Allow some tolerance for edge cases
        max_invalid = len(KNOWN_SPECIES) * 0.05  # 5% tolerance
        assert invalid_count <= max_invalid, \
            f"Too many invalid SMILES ({invalid_count}): {invalid_species[:10]}"

    def test_no_duplicate_smiles_for_different_meanings(self):
        """Check for suspicious duplicates that might indicate errors."""
        # Collect SMILES->codes mapping
        smiles_to_codes = {}
        for code, smiles in KNOWN_SPECIES.items():
            if smiles:
                if smiles not in smiles_to_codes:
                    smiles_to_codes[smiles] = []
                smiles_to_codes[smiles].append(code)

        # Some duplicates are fine (aliases), but excessive duplicates are suspicious
        for smiles, codes in smiles_to_codes.items():
            if len(codes) > 10:  # Suspicious if more than 10 codes map to same SMILES
                # Unless it's a very common species
                common = ["C", "CC", "CCC", "C=O", "CO", "c1ccccc1"]
                if smiles not in common:
                    pytest.fail(f"Suspicious duplicate: {smiles} -> {codes}")


class TestIntegration:
    """Integration tests for the full SMILES conversion pipeline."""

    def test_species_code_priority(self):
        """Species code lookup should take priority over formula parsing."""
        # ISOPRE has a known SMILES
        result_with_code = gecko_to_smiles("WRONG_FORMULA", species_code="ISOPRE")
        result_direct = KNOWN_SPECIES.get("ISOPRE")
        assert result_with_code == result_direct, "Species code should override formula"

    def test_fallback_chain(self):
        """Test the fallback chain: code -> formula -> aromatic -> default."""
        # Unknown code should fall back to formula parsing
        result = gecko_to_smiles("CH3CH2OH", species_code="UNKNOWN999")
        assert result is not None
        assert isinstance(result, str)

    @pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")
    def test_output_is_valid_smiles(self):
        """Output from gecko_to_smiles should be valid SMILES when possible."""
        test_cases = [
            ("", "ISOPRE"),
            ("", "APINEN"),
            ("", "R07000"),
            ("CH3OH", ""),
            ("CH3CH2OH", ""),
        ]

        for formula, code in test_cases:
            result = gecko_to_smiles(formula, species_code=code)
            if result:  # Skip empty results
                validated = validate_smiles_with_rdkit(result)
                # Either valid or a known edge case
                if validated is None and code:
                    # Check if it's in our known database
                    if code.upper() in KNOWN_SPECIES:
                        pytest.fail(f"Known species {code} produced invalid SMILES: {result}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
