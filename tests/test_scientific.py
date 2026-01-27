"""
Scientific validation tests for GECKO-A Web Interface.

These tests verify the correctness of core scientific calculations:
1. C* (saturation concentration) calculation
2. Gas-particle partitioning (Fp) calculation
3. Vapor pressure reading and conversion
4. Unit conversions
5. Iterative partitioning convergence

References:
- Pankow (1994): Atmos. Environ. 28, 185-188
- Donahue et al. (2006): Environ. Sci. Technol. 40, 2635-2643
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from gecko_web.postprocessing import (
    PartitioningCalculator,
    VaporPressureReader,
    FortranOutputParser,
    _estimate_mw_from_composition,
)


class TestCStarCalculation:
    """Test saturation concentration (C*) calculations."""

    def test_cstar_formula_validation(self):
        """
        Verify C* calculation against hand-calculated values.

        C* = (MW * Pvap_Pa) / (R * T) * 10^6  [µg/m³]

        For MW = 100 g/mol, Pvap = 1e-6 atm (log = -6), T = 298 K:
        Pvap_Pa = 1e-6 * 101325 = 0.101325 Pa
        C* = (100 * 0.101325) / (8.314 * 298) * 1e6
           = 10.1325 / 2477.572 * 1e6
           = 0.00409 * 1e6
           = 4090 µg/m³
        """
        calc = PartitioningCalculator(temperature_k=298.0)
        c_star = calc.calculate_c_star(mw=100, log_pvap_atm=-6)

        # Hand calculation
        expected = (100 * 1e-6 * 101325) / (8.314 * 298) * 1e6

        assert abs(c_star - expected) / expected < 0.01, \
            f"C* mismatch: got {c_star:.2f}, expected {expected:.2f}"

    def test_cstar_temperature_dependence(self):
        """Verify C* increases with temperature (inverse relationship via ideal gas)."""
        calc_cold = PartitioningCalculator(temperature_k=273.0)
        calc_hot = PartitioningCalculator(temperature_k=313.0)

        c_star_cold = calc_cold.calculate_c_star(mw=100, log_pvap_atm=-6)
        c_star_hot = calc_hot.calculate_c_star(mw=100, log_pvap_atm=-6)

        # At higher T, denominator (R*T) increases, so C* should decrease
        # But Pvap typically increases faster with T - this test uses fixed Pvap
        # so C* should decrease with T
        assert c_star_cold > c_star_hot, \
            "C* should decrease with temperature at fixed Pvap"

    def test_cstar_volatility_range(self):
        """
        Test C* spans reasonable range for atmospheric species.

        VOC: log(C*) > 6 (very volatile)
        IVOC: log(C*) 3-6
        SVOC: log(C*) 0-3
        LVOC: log(C*) -3 to 0
        ELVOC: log(C*) < -3
        """
        calc = PartitioningCalculator(temperature_k=298.0)

        # VOC (high Pvap)
        c_star_voc = calc.calculate_c_star(mw=100, log_pvap_atm=-2)
        assert np.log10(c_star_voc) > 4, "VOC should have log(C*) > 4"

        # SVOC (medium Pvap)
        c_star_svoc = calc.calculate_c_star(mw=200, log_pvap_atm=-8)
        assert 0 < np.log10(c_star_svoc) < 4, "SVOC should have log(C*) 0-4"

        # LVOC (low Pvap)
        c_star_lvoc = calc.calculate_c_star(mw=300, log_pvap_atm=-12)
        assert np.log10(c_star_lvoc) < 1, "LVOC should have log(C*) < 1"

    def test_cstar_edge_cases(self):
        """Test edge cases in C* calculation."""
        calc = PartitioningCalculator()

        # Invalid MW should return very high C* (very volatile)
        c_star = calc.calculate_c_star(mw=0, log_pvap_atm=-6)
        assert c_star == 1e10, "Zero MW should return 1e10"

        # Very high Pvap (log > 10) should return very high C*
        c_star = calc.calculate_c_star(mw=100, log_pvap_atm=11)
        assert c_star == 1e10, "Very high Pvap should return 1e10"


class TestPartitioningCalculation:
    """Test gas-particle partitioning (Fp) calculations."""

    def test_fp_formula_validation(self):
        """
        Verify Fp = C_OA / (C_OA + C*).

        For C_OA = 10 µg/m³, C* = 10 µg/m³:
        Fp = 10 / (10 + 10) = 0.5
        """
        calc = PartitioningCalculator(seed_mass=10.0)

        species_data = [{
            'code': 'TEST',
            'mw': 100,
            'c_star': 10.0,  # C* = 10 µg/m³
            'concentration': 1.0,  # Low concentration won't affect C_OA much
        }]

        result = calc.calculate_partitioning_static(species_data, c_oa_assumed=10.0)

        expected_fp = 10.0 / (10.0 + 10.0)  # 0.5
        assert abs(result.iloc[0]['f_p'] - expected_fp) < 0.01, \
            f"Fp mismatch: got {result.iloc[0]['f_p']:.3f}, expected {expected_fp:.3f}"

    def test_fp_limits(self):
        """
        Test Fp approaches limits:
        - Fp → 1 when C* << C_OA (low volatility, mostly in particle)
        - Fp → 0 when C* >> C_OA (high volatility, mostly in gas)
        """
        calc = PartitioningCalculator(seed_mass=10.0)

        # Low volatility (C* = 0.001, C_OA = 10)
        species_low_vol = [{'code': 'LOW', 'mw': 100, 'c_star': 0.001, 'concentration': 0.1}]
        result_low = calc.calculate_partitioning_static(species_low_vol, c_oa_assumed=10.0)
        assert result_low.iloc[0]['f_p'] > 0.999, \
            f"Low volatility species should have Fp ≈ 1, got {result_low.iloc[0]['f_p']:.4f}"

        # High volatility (C* = 100000, C_OA = 10)
        species_high_vol = [{'code': 'HIGH', 'mw': 100, 'c_star': 100000, 'concentration': 0.1}]
        result_high = calc.calculate_partitioning_static(species_high_vol, c_oa_assumed=10.0)
        assert result_high.iloc[0]['f_p'] < 0.001, \
            f"High volatility species should have Fp ≈ 0, got {result_high.iloc[0]['f_p']:.4f}"

    def test_dynamic_partitioning_convergence(self):
        """Test that iterative partitioning converges."""
        calc = PartitioningCalculator(seed_mass=1.0, temperature_k=298.0)

        # Create test species with range of volatilities
        species_data = [
            {'code': 'SP1', 'mw': 150, 'c_star': 0.1, 'concentration': 1.0},
            {'code': 'SP2', 'mw': 200, 'c_star': 1.0, 'concentration': 2.0},
            {'code': 'SP3', 'mw': 250, 'c_star': 10.0, 'concentration': 3.0},
            {'code': 'SP4', 'mw': 300, 'c_star': 100.0, 'concentration': 5.0},
        ]

        result = calc.calculate_partitioning_dynamic(species_data, max_iterations=100)

        # Should not return empty DataFrame
        assert len(result) == 4, "Should return results for all 4 species"

        # Fp should be between 0 and 1
        assert all(0 <= fp <= 1 for fp in result['f_p']), "All Fp values should be in [0, 1]"

        # Particle + gas should equal total
        for _, row in result.iterrows():
            total = row['c_particle'] + row['c_gas']
            assert abs(total - row['c_total']) < 1e-6, \
                f"Mass balance: particle + gas should equal total"

    def test_dynamic_vs_static_with_seed(self):
        """
        With high seed mass, dynamic should approach static solution.

        When seed >> organic contribution, C_OA ≈ seed, and dynamic becomes static.
        """
        seed_mass = 100.0
        calc = PartitioningCalculator(seed_mass=seed_mass)

        # Low organic concentrations relative to seed
        species_data = [
            {'code': 'SP1', 'mw': 200, 'c_star': 10.0, 'concentration': 0.01},
        ]

        result_dynamic = calc.calculate_partitioning_dynamic(species_data)
        result_static = calc.calculate_partitioning_static(species_data, c_oa_assumed=seed_mass)

        # Fp should be very similar
        fp_dynamic = result_dynamic.iloc[0]['f_p']
        fp_static = result_static.iloc[0]['f_p']

        assert abs(fp_dynamic - fp_static) < 0.1, \
            f"With large seed, dynamic ({fp_dynamic:.3f}) should ≈ static ({fp_static:.3f})"

    def test_mass_conservation(self):
        """Verify mass is conserved in partitioning."""
        calc = PartitioningCalculator(seed_mass=5.0)

        species_data = [
            {'code': 'A', 'mw': 150, 'c_star': 1.0, 'concentration': 10.0},
            {'code': 'B', 'mw': 200, 'c_star': 10.0, 'concentration': 20.0},
            {'code': 'C', 'mw': 250, 'c_star': 100.0, 'concentration': 15.0},
        ]

        result = calc.calculate_partitioning_dynamic(species_data)

        total_input = sum(sp['concentration'] for sp in species_data)
        total_output = result['c_particle'].sum() + result['c_gas'].sum()

        assert abs(total_input - total_output) < 1e-6, \
            f"Mass not conserved: input={total_input:.4f}, output={total_output:.4f}"


class TestVaporPressureEstimation:
    """Test vapor pressure estimation fallback."""

    def test_mw_correlation(self):
        """Test MW-based Pvap estimation."""
        # Higher MW should give lower Pvap
        log_pvap_low_mw = VaporPressureReader.estimate_from_mw(100)
        log_pvap_high_mw = VaporPressureReader.estimate_from_mw(300)

        assert log_pvap_low_mw > log_pvap_high_mw, \
            "Higher MW should have lower vapor pressure"

    def test_oxygen_reduces_volatility(self):
        """Oxygen atoms should reduce volatility (lower Pvap)."""
        log_pvap_no_o = VaporPressureReader.estimate_from_mw(150, n_oxygen=0)
        log_pvap_with_o = VaporPressureReader.estimate_from_mw(150, n_oxygen=4)

        assert log_pvap_no_o > log_pvap_with_o, \
            "Oxygenated species should have lower Pvap"

        # Each O should reduce by ~0.8 orders of magnitude
        expected_reduction = 4 * 0.8
        actual_reduction = log_pvap_no_o - log_pvap_with_o
        assert abs(actual_reduction - expected_reduction) < 0.1, \
            f"Oxygen effect: expected {expected_reduction:.1f}, got {actual_reduction:.1f}"

    def test_invalid_mw(self):
        """Invalid MW should return 0."""
        log_pvap = VaporPressureReader.estimate_from_mw(0)
        assert log_pvap == 0.0, "Zero MW should return log(Pvap) = 0"


class TestFortranParser:
    """Test Fortran output parsing."""

    def test_smart_split_merged_columns(self):
        """Test splitting merged scientific notation columns."""
        # Typical Fortran output where columns run together
        line = "1.00E+012.00E+013.50E-05"
        parts = FortranOutputParser.smart_split(line)

        assert len(parts) == 3, f"Should split into 3 parts, got {len(parts)}"
        assert float(parts[0]) == pytest.approx(1.0e1)
        assert float(parts[1]) == pytest.approx(2.0e1)
        assert float(parts[2]) == pytest.approx(3.5e-5)

    def test_smart_split_spaced_columns(self):
        """Test normal space-separated columns."""
        line = "SPECIES  1.23E+05  4.56E-03  789.0"
        parts = FortranOutputParser.smart_split(line)

        assert len(parts) == 4
        assert parts[0] == "SPECIES"
        assert float(parts[1]) == pytest.approx(1.23e5)

    def test_smart_split_negative_exponents(self):
        """Test handling of negative exponents."""
        line = "-1.23E-05-4.56E+02"
        parts = FortranOutputParser.smart_split(line)

        # Should handle negative numbers correctly
        assert len(parts) == 2
        assert float(parts[0]) == pytest.approx(-1.23e-5)
        assert float(parts[1]) == pytest.approx(-4.56e2)


class TestMolecularWeightEstimation:
    """Test molecular weight estimation from composition."""

    def test_mw_from_composition(self):
        """Test MW calculation from element counts."""
        # Methane: CH4
        # C=12, H=4 -> MW = 16
        mw = _estimate_mw_from_composition(n_c=1, n_h=4, n_n=0, n_o=0)
        assert abs(mw - 16.0) < 0.1, f"CH4 MW should be 16, got {mw}"

        # Water: H2O
        # H=2, O=16 -> MW = 18
        mw = _estimate_mw_from_composition(n_c=0, n_h=2, n_n=0, n_o=1)
        assert abs(mw - 18.0) < 0.1, f"H2O MW should be 18, got {mw}"

        # Nitric acid: HNO3
        # H=1, N=14, O=48 -> MW = 63
        mw = _estimate_mw_from_composition(n_c=0, n_h=1, n_n=1, n_o=3)
        assert abs(mw - 63.0) < 0.1, f"HNO3 MW should be 63, got {mw}"


class TestScientificRealism:
    """Test that outputs are scientifically realistic."""

    def test_soa_yield_reasonable(self):
        """SOA yield should be between 0 and 1."""
        calc = PartitioningCalculator(seed_mass=10.0)

        # Simulate oxidation products
        species_data = [
            {'code': f'P{i}', 'mw': 150 + i*20, 'c_star': 10**(i-2), 'concentration': 1.0}
            for i in range(10)
        ]

        result = calc.calculate_partitioning_dynamic(species_data)

        total_organic = result['c_total'].sum()
        particle_organic = result['c_particle'].sum()

        if total_organic > 0:
            yield_value = particle_organic / total_organic
            assert 0 <= yield_value <= 1, \
                f"SOA yield should be in [0, 1], got {yield_value:.3f}"

    def test_vbs_distribution_shape(self):
        """VBS distribution should show expected pattern."""
        calc = PartitioningCalculator(seed_mass=5.0)

        # Create species spanning VBS bins from -3 to 6
        species_data = []
        for log_cstar in range(-3, 7):
            species_data.append({
                'code': f'BIN{log_cstar}',
                'mw': 200,
                'c_star': 10**log_cstar,
                'concentration': 10.0,  # Same total for each bin
            })

        result = calc.calculate_partitioning_dynamic(species_data)
        result['log_c_star'] = np.log10(result['c_star'])

        # Low volatility bins should have higher Fp
        low_vol = result[result['log_c_star'] < 0]['f_p'].mean()
        high_vol = result[result['log_c_star'] > 3]['f_p'].mean()

        assert low_vol > high_vol, \
            f"Low volatility species should have higher Fp: low={low_vol:.3f}, high={high_vol:.3f}"


class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_partitioning_workflow(self):
        """Test complete partitioning workflow from Pvap to Fp."""
        calc = PartitioningCalculator(temperature_k=298.0, seed_mass=5.0)

        # Create realistic species with estimated vapor pressures
        test_species = []
        for i, (mw, n_o) in enumerate([(100, 1), (150, 2), (200, 3), (250, 4)]):
            log_pvap = VaporPressureReader.estimate_from_mw(mw, n_oxygen=n_o)
            c_star = calc.calculate_c_star(mw, log_pvap)

            test_species.append({
                'code': f'TEST{i}',
                'mw': mw,
                'c_star': c_star,
                'concentration': 1.0 + i * 0.5,
            })

        result = calc.calculate_partitioning_dynamic(test_species)

        # Verify results are reasonable
        assert len(result) == 4
        assert all(result['f_p'] >= 0)
        assert all(result['f_p'] <= 1)
        assert all(result['c_star'] > 0)

        # Higher MW + more oxygen should give lower C* and higher Fp
        fp_values = result['f_p'].values
        # Generally should increase with index (lower volatility)
        # Allow some variability
        assert fp_values[-1] > fp_values[0], \
            "Higher oxygenation should lead to higher Fp"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
