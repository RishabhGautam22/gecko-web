import numpy as np
import pandas as pd

from gecko_web.postprocessing import PartitioningCalculator


def test_c_star_formula():
    calc = PartitioningCalculator(temperature_k=298.0)
    mw = 100.0  # g/mol
    log_pvap_atm = -6.0  # 1e-6 atm
    pvap_pa = 10 ** log_pvap_atm * 101325.0
    expected = (mw * pvap_pa) / (calc.GAS_CONSTANT * 298.0) * 1e6

    c_star = calc.calculate_c_star(mw, log_pvap_atm)
    assert np.isclose(c_star, expected, rtol=1e-6)


def test_partitioning_converges_to_closed_form():
    # Single species: C* = 10, C_total = 100, seed = 0
    # Closed-form solution: C_OA = C_total - C* = 90
    calc = PartitioningCalculator(temperature_k=298.0, seed_mass=0.0)
    df = calc.calculate_partitioning_dynamic([
        {'code': 'X', 'mw': 100.0, 'c_star': 10.0, 'concentration': 100.0}
    ])
    c_oa = float(df['c_oa_equilibrium'].iloc[0])
    f_p = float(df['f_p'].iloc[0])

    assert np.isclose(c_oa, 90.0, rtol=1e-3)
    assert np.isclose(f_p, 0.9, rtol=1e-3)


def test_seed_mass_sets_minimum_c_oa():
    # Very volatile species with low partitioning should still respect seed mass
    calc = PartitioningCalculator(temperature_k=298.0, seed_mass=20.0)
    df = calc.calculate_partitioning_dynamic([
        {'code': 'X', 'mw': 100.0, 'c_star': 1e6, 'concentration': 100.0}
    ])
    c_oa = float(df['c_oa_equilibrium'].iloc[0])
    assert c_oa >= 20.0
