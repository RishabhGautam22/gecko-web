from gecko_web.main import _reduce_pfa


def test_pfa_keeps_precursors_for_target():
    species = ['A', 'B', 'C']
    reactions = [
        {'reactants': ['A'], 'products': ['B']},
        {'reactants': ['B'], 'products': ['C']},
    ]

    kept_species, _ = _reduce_pfa(
        species=species,
        reactions=reactions,
        target_species=['C'],
        error_threshold=0.01,
        min_species=0
    )

    assert 'A' in kept_species
    assert 'B' in kept_species
    assert 'C' in kept_species


def test_pfa_excludes_downstream_for_upstream_target():
    species = ['A', 'B', 'C']
    reactions = [
        {'reactants': ['A'], 'products': ['B']},
        {'reactants': ['B'], 'products': ['C']},
    ]

    kept_species, _ = _reduce_pfa(
        species=species,
        reactions=reactions,
        target_species=['A'],
        error_threshold=0.01,
        min_species=0
    )

    assert kept_species == ['A'] or set(kept_species) == {'A'}
