"""
GECKO-A Chemical Database Package

This package provides comprehensive chemical data for atmospheric chemistry modeling,
including verified SMILES, molecular properties, and reaction parameters.

All data has been validated against:
- NIST Chemistry WebBook
- PubChem
- EPA CompTox Dashboard
- IUPAC Gold Book
- Primary literature sources

Author: GECKO-A Development Team
"""

from .compound_database import (
    CompoundDatabase,
    Compound,
    COMPOUND_DATABASE,
    get_database,
    get_compound,
    get_all_compounds,
    get_compounds_by_category,
    get_gecko_formula,
    get_smiles,
    validate_compound,
    search_compounds,
    get_atmospheric_lifetime,
    get_compound_names_for_dropdown
)

from .voc_categories import (
    VOC_CATEGORIES,
    VOCCategory,
    VOCSubcategory,
    get_category_compounds,
    get_subcategory_compounds,
    get_all_compounds as get_all_voc_compounds,
    get_high_soa_compounds
)

from .reaction_data import (
    ReactionDatabase,
    ReactionKinetics,
    ArrheniusParams,
    TroePressureParams,
    OxidantType,
    ReactionType,
    reaction_database,
    get_rate_constant,
    get_branching_ratios,
    estimate_soa_yield,
    calculate_atmospheric_lifetime
)

__all__ = [
    # Compound database
    'CompoundDatabase',
    'Compound',
    'COMPOUND_DATABASE',
    'get_database',
    'get_compound',
    'get_all_compounds',
    'get_compounds_by_category',
    'get_gecko_formula',
    'get_smiles',
    'validate_compound',
    'search_compounds',
    'get_atmospheric_lifetime',
    'get_compound_names_for_dropdown',
    # VOC categories
    'VOC_CATEGORIES',
    'VOCCategory',
    'VOCSubcategory',
    'get_category_compounds',
    'get_subcategory_compounds',
    'get_all_voc_compounds',
    'get_high_soa_compounds',
    # Reaction kinetics
    'ReactionDatabase',
    'ReactionKinetics',
    'ArrheniusParams',
    'TroePressureParams',
    'OxidantType',
    'ReactionType',
    'reaction_database',
    'get_rate_constant',
    'get_branching_ratios',
    'estimate_soa_yield',
    'calculate_atmospheric_lifetime'
]
