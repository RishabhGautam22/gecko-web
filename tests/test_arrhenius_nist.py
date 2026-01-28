"""
Arrhenius Rate Constant Validation Tests Against NIST/JPL/IUPAC Reference Data

This test module validates the implementation of Arrhenius rate constant calculations
against experimentally determined values from authoritative sources:

- NIST Chemical Kinetics Database
- NASA/JPL Publication 19-5 (2020): Chemical Kinetics and Photochemical Data
- IUPAC Task Group on Atmospheric Chemical Kinetic Data Evaluation

Test Philosophy:
- Each test compares calculated k values against literature values at specific temperatures
- Acceptance criterion: within 10% of reference value (typical experimental uncertainty)
- Tests cover both simple Arrhenius and Troe pressure-dependent reactions

References:
- Burkholder, J.B., et al. (2020) JPL Publication 19-5
- Atkinson, R., et al. (2006) ACP 6, 3625-4055
- Sander, S.P., et al. (2011) JPL Publication 10-6

Author: GECKO-A Development Team
"""

import pytest
import math
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from gecko_web.chemdata.reaction_data import (
    ArrheniusParams,
    TroePressureParams,
    R_GAS,
)


# ==============================================================================
# NIST/JPL Reference Data
# ==============================================================================

# Structure: (reaction_name, {params}, {reference_k_at_T})
NIST_REFERENCE_DATA = {
    # -------------------------------------------------------------------------
    # OH + VOC Reactions (JPL 19-5, Table 1-1)
    # -------------------------------------------------------------------------
    "OH + CH4": {
        "arrhenius": ArrheniusParams(
            A=2.45e-12,
            n=0,
            Ea=14818.0,  # J/mol (1780 K * R)
            reference="JPL 19-5 (2020)"
        ),
        "reference_values": {
            # Temperature (K): k (cm³ molecule⁻¹ s⁻¹) calculated from JPL parameters
            # Using k = A * exp(-Ea/RT) with A=2.45e-12, Ea/R=1780K
            298: 6.31e-15,  # JPL recommended value at 298K
        },
        "tolerance": 0.10,  # 10% tolerance for self-consistency check
        "source": "JPL 19-5 Table 1-1; self-consistency test"
    },

    "OH + C2H6": {
        "arrhenius": ArrheniusParams(
            A=7.66e-12,
            n=0,
            Ea=8650.0,  # 1040 K * R
            reference="JPL 19-5 (2020)"
        ),
        "reference_values": {
            250: 1.13e-13,
            298: 2.41e-13,
            350: 4.42e-13,
        },
        "tolerance": 0.15,
        "source": "JPL 19-5 Table 1-1"
    },

    "OH + C3H8": {
        "arrhenius": ArrheniusParams(
            A=7.60e-12,
            n=0,
            Ea=5820.0,  # 700 K * R
            reference="JPL 19-5 (2020)"
        ),
        "reference_values": {
            # Calculated from JPL parameters: k = 7.6e-12 * exp(-Ea/(R*T))
            # where Ea = 5820 J/mol, R = 8.314 J/(mol*K), T = 298 K
            # k = 7.6e-12 * exp(-5820/(8.314*298)) = 7.6e-12 * exp(-2.35) = 7.25e-13
            298: 7.25e-13,  # Self-consistency check
        },
        "tolerance": 0.05,  # Tight tolerance for self-consistency
        "source": "JPL 19-5 Table 1-1; self-consistency test"
    },

    "OH + CO": {
        "arrhenius": ArrheniusParams(
            A=1.44e-13,
            n=0,
            Ea=0.0,
            reference="JPL 19-5 (2020), low-pressure limit"
        ),
        "reference_values": {
            # Note: This is the low-pressure limit; actual k is pressure-dependent
            298: 1.44e-13,
        },
        "tolerance": 0.20,  # Higher tolerance due to pressure dependence
        "source": "JPL 19-5 (k at low pressure)"
    },

    # -------------------------------------------------------------------------
    # O3 + Alkene Reactions (IUPAC/Atkinson 2006)
    # -------------------------------------------------------------------------
    "O3 + C2H4": {
        "arrhenius": ArrheniusParams(
            A=9.1e-15,
            n=0,
            Ea=21560.0,  # 2594 K * R
            reference="IUPAC (2006)"
        ),
        "reference_values": {
            # Calculated from IUPAC parameters: k = 9.1e-15 * exp(-2594/T)
            298: 1.59e-18,  # Self-consistency check at 298K
        },
        "tolerance": 0.10,
        "source": "IUPAC Sheet Ox_VOC6; self-consistency test"
    },

    "O3 + C3H6": {
        "arrhenius": ArrheniusParams(
            A=5.5e-15,
            n=0,
            Ea=15750.0,  # 1894 K * R
            reference="IUPAC (2006)"
        ),
        "reference_values": {
            # Calculated from IUPAC parameters: k = 5.5e-15 * exp(-1894/T)
            298: 1.01e-17,  # Self-consistency at 298K
        },
        "tolerance": 0.10,
        "source": "IUPAC Sheet Ox_VOC8; self-consistency test"
    },

    # -------------------------------------------------------------------------
    # Isoprene Chemistry (MCM v3.3.1, Carlton 2009)
    # -------------------------------------------------------------------------
    "OH + ISOPRENE": {
        "arrhenius": ArrheniusParams(
            A=2.54e-11,
            n=0,
            Ea=-3300.0,  # Negative Ea (rate increases with T)
            reference="MCM v3.3.1"
        ),
        "reference_values": {
            298: 1.0e-10,  # ~1.0e-10 at 298 K
        },
        "tolerance": 0.15,
        "source": "MCM v3.3.1; Atkinson & Arey (2003)"
    },

    # -------------------------------------------------------------------------
    # NO3 Reactions (JPL 19-5)
    # -------------------------------------------------------------------------
    "NO3 + NO2": {
        "arrhenius": ArrheniusParams(
            A=1.0e-12,
            n=0,
            Ea=0.0,
            reference="JPL 19-5 (equilibrium)"
        ),
        "reference_values": {
            298: 1.0e-12,
        },
        "tolerance": 0.30,  # Higher uncertainty for this reaction
        "source": "JPL 19-5 (part of N2O5 equilibrium)"
    },
}


# ==============================================================================
# Troe Pressure-Dependent Reactions Reference Data
# ==============================================================================

TROE_REFERENCE_DATA = {
    # -------------------------------------------------------------------------
    # OH + NO2 + M -> HNO3 + M (JPL 19-5)
    # -------------------------------------------------------------------------
    "OH + NO2 -> HNO3": {
        "troe": TroePressureParams(
            k0_A=1.8e-30,
            k0_n=-3.0,
            k0_Ea=0,
            kinf_A=2.8e-11,
            kinf_n=0,
            kinf_Ea=0,
            Fc=0.6,
            reference="JPL 19-5"
        ),
        "reference_values": {
            # (T, M): k
            (298, 2.5e19): 1.1e-11,  # ~1 atm
            (250, 2.5e19): 1.3e-11,
        },
        "tolerance": 0.25,
        "source": "JPL 19-5 Table 2-1"
    },

    # -------------------------------------------------------------------------
    # HO2 + NO2 + M -> HNO4 + M (JPL 19-5)
    # -------------------------------------------------------------------------
    "HO2 + NO2 -> HNO4": {
        "troe": TroePressureParams(
            k0_A=1.9e-31,
            k0_n=-3.4,
            k0_Ea=0,
            kinf_A=4.0e-12,
            kinf_n=-0.3,
            kinf_Ea=0,
            Fc=0.6,
            reference="JPL 19-5"
        ),
        "reference_values": {
            (298, 2.5e19): 1.4e-12,
            (250, 2.5e19): 2.5e-12,
        },
        "tolerance": 0.30,
        "source": "JPL 19-5 Table 2-1"
    },
}


# ==============================================================================
# Test Classes
# ==============================================================================

class TestArrheniusCalculation:
    """Test Arrhenius rate constant calculation against NIST reference values."""

    @pytest.mark.parametrize("reaction_name,data", NIST_REFERENCE_DATA.items())
    def test_arrhenius_vs_nist(self, reaction_name, data):
        """
        Verify Arrhenius k calculation matches NIST/JPL reference values.

        Test criterion: |k_calc - k_ref| / k_ref < tolerance
        """
        params = data["arrhenius"]
        tolerance = data["tolerance"]
        source = data["source"]

        for T, k_ref in data["reference_values"].items():
            k_calc = params.calculate_k(T)

            rel_error = abs(k_calc - k_ref) / k_ref

            assert rel_error < tolerance, (
                f"Arrhenius validation FAILED for {reaction_name} at T={T} K\n"
                f"  Calculated: k = {k_calc:.2e} cm³/molec/s\n"
                f"  Reference:  k = {k_ref:.2e} cm³/molec/s\n"
                f"  Relative error: {rel_error*100:.1f}% (tolerance: {tolerance*100:.0f}%)\n"
                f"  Source: {source}"
            )

    def test_temperature_dependence_direction(self):
        """
        Verify physical consistency: positive Ea → k decreases with T decrease.
        """
        # CH4 + OH: positive Ea, endothermic barrier
        params = NIST_REFERENCE_DATA["OH + CH4"]["arrhenius"]

        k_298 = params.calculate_k(298)
        k_250 = params.calculate_k(250)
        k_350 = params.calculate_k(350)

        assert k_250 < k_298 < k_350, (
            f"Temperature dependence wrong for positive Ea reaction:\n"
            f"  k(250K) = {k_250:.2e} should be < k(298K) = {k_298:.2e} < k(350K) = {k_350:.2e}"
        )

    def test_negative_ea_increases_with_temperature(self):
        """
        Verify physical consistency: negative Ea → k increases with T decrease.
        """
        # ISOPRENE + OH: negative Ea (addition reaction, no barrier)
        params = NIST_REFERENCE_DATA["OH + ISOPRENE"]["arrhenius"]

        k_298 = params.calculate_k(298)
        k_250 = params.calculate_k(250)

        # For negative Ea, k should increase at lower T
        assert k_250 > k_298, (
            f"Temperature dependence wrong for negative Ea reaction:\n"
            f"  k(250K) = {k_250:.2e} should be > k(298K) = {k_298:.2e}"
        )

    def test_uncertainty_bounds(self):
        """Verify uncertainty bounds bracket the reference value."""
        for reaction_name, data in NIST_REFERENCE_DATA.items():
            params = data["arrhenius"]

            for T, k_ref in data["reference_values"].items():
                k_low, k_high = params.calculate_k_uncertainty(T)

                # Reference value should be within uncertainty bounds
                # (with some margin for experimental uncertainty in reference)
                margin = 1.5  # Allow 50% margin for experimental scatter

                assert k_low / margin <= k_ref <= k_high * margin, (
                    f"Reference value outside uncertainty bounds for {reaction_name} at T={T}K\n"
                    f"  Reference: {k_ref:.2e}\n"
                    f"  Bounds: [{k_low:.2e}, {k_high:.2e}]"
                )


class TestTroeCalculation:
    """Test Troe pressure-dependent rate constant calculations."""

    @pytest.mark.parametrize("reaction_name,data", TROE_REFERENCE_DATA.items())
    def test_troe_vs_reference(self, reaction_name, data):
        """
        Verify Troe k calculation matches JPL reference values.
        """
        params = data["troe"]
        tolerance = data["tolerance"]
        source = data["source"]

        for (T, M), k_ref in data["reference_values"].items():
            k_calc = params.calculate_k(T, M)

            rel_error = abs(k_calc - k_ref) / k_ref

            assert rel_error < tolerance, (
                f"Troe validation FAILED for {reaction_name} at T={T}K, M={M:.1e}\n"
                f"  Calculated: k = {k_calc:.2e}\n"
                f"  Reference:  k = {k_ref:.2e}\n"
                f"  Relative error: {rel_error*100:.1f}% (tolerance: {tolerance*100:.0f}%)\n"
                f"  Source: {source}"
            )

    def test_pressure_dependence(self):
        """
        Verify Troe falloff behavior:
        - At low pressure: k ∝ [M]
        - At high pressure: k → k∞
        """
        params = TROE_REFERENCE_DATA["OH + NO2 -> HNO3"]["troe"]
        T = 298

        M_low = 1e17  # Low pressure
        M_mid = 2.5e19  # ~1 atm
        M_high = 1e22  # High pressure

        k_low = params.calculate_k(T, M_low)
        k_mid = params.calculate_k(T, M_mid)
        k_high = params.calculate_k(T, M_high)

        # At low pressure, k should scale roughly with M
        # At high pressure, k should approach k_inf
        assert k_low < k_mid < k_high, (
            f"Pressure dependence incorrect:\n"
            f"  k(low P) = {k_low:.2e}\n"
            f"  k(1 atm) = {k_mid:.2e}\n"
            f"  k(high P) = {k_high:.2e}"
        )

        # High pressure limit check
        k_inf = params.kinf_A * (T / 300.0) ** params.kinf_n
        assert k_high / k_inf > 0.8, (
            f"k at high pressure should approach k_inf:\n"
            f"  k(high P) = {k_high:.2e}\n"
            f"  k_inf = {k_inf:.2e}\n"
            f"  Ratio: {k_high/k_inf:.2f}"
        )


class TestPhysicalConstants:
    """Verify physical constants match standard values."""

    def test_gas_constant(self):
        """R should be 8.314 J/(mol·K) (CODATA 2018)."""
        assert abs(R_GAS - 8.314) < 0.001, f"R = {R_GAS}, expected 8.314"

    def test_arrhenius_at_300K(self):
        """At T=300K, (T/300)^n = 1 for any n."""
        params = ArrheniusParams(A=1e-12, n=2.5, Ea=0)
        k = params.calculate_k(300)
        assert abs(k - 1e-12) < 1e-20, "At 300K, k should equal A when Ea=0"

    def test_zero_ea_gives_constant_k(self):
        """With Ea=0 and n=0, k should be constant with T."""
        params = ArrheniusParams(A=5e-12, n=0, Ea=0)

        k_250 = params.calculate_k(250)
        k_300 = params.calculate_k(300)
        k_350 = params.calculate_k(350)

        assert k_250 == k_300 == k_350 == 5e-12, (
            f"k should be constant when Ea=0 and n=0:\n"
            f"  k(250K) = {k_250:.2e}\n"
            f"  k(300K) = {k_300:.2e}\n"
            f"  k(350K) = {k_350:.2e}"
        )


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_low_temperature(self):
        """Rate constant should still be positive at low T."""
        params = ArrheniusParams(A=1e-12, n=0, Ea=10000)

        k = params.calculate_k(150)  # Very cold stratosphere
        assert k > 0, f"k must be positive: {k}"
        assert math.isfinite(k), f"k must be finite: {k}"

    def test_very_high_temperature(self):
        """Rate constant should be reasonable at high T."""
        params = ArrheniusParams(A=1e-12, n=0, Ea=10000)

        k = params.calculate_k(500)  # Hot conditions
        assert k > 0, f"k must be positive: {k}"
        assert k < 1e-6, f"k should not be unreasonably large: {k}"

    def test_large_negative_ea(self):
        """Large negative Ea should give reasonable k at low T."""
        params = ArrheniusParams(A=1e-12, n=0, Ea=-5000)  # Strong T⁻ dependence

        k_200 = params.calculate_k(200)
        k_300 = params.calculate_k(300)

        assert k_200 > k_300, "Negative Ea: k should increase at lower T"
        assert k_200 < 1e-6, f"k should not be unreasonably large: {k_200}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
