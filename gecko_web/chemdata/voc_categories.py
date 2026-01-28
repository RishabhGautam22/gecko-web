"""
VOC Categories and Classification

This module provides hierarchical categorization of volatile organic compounds (VOCs)
for atmospheric chemistry modeling with GECKO-A.

Categories are based on:
- Chemical structure and functional groups
- Emission sources (biogenic, anthropogenic)
- Atmospheric reactivity
- SOA formation potential

References:
- Goldstein & Galbally (2007): Known and unexplored organic constituents in the atmosphere
- Atkinson & Arey (2003): Atmospheric degradation of volatile organic compounds
- Guenther et al. (2012): MEGAN biogenic emissions

Author: GECKO-A Development Team
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


# ==============================================================================
# VOC Category Definitions
# ==============================================================================

@dataclass
class VOCSubcategory:
    """Represents a subcategory of VOCs."""
    name: str
    description: str
    compounds: List[str]  # List of compound names
    typical_sources: List[str]
    soa_potential: str  # 'low', 'medium', 'high', 'very_high'


@dataclass
class VOCCategory:
    """Represents a main category of VOCs."""
    name: str
    description: str
    subcategories: Dict[str, VOCSubcategory]
    emission_type: str  # 'biogenic', 'anthropogenic', 'both'


# ==============================================================================
# Complete VOC Categories Database
# ==============================================================================

VOC_CATEGORIES = {
    # =========================================================================
    # ALKANES (Saturated Hydrocarbons)
    # =========================================================================
    'alkanes': VOCCategory(
        name='Alkanes',
        description='Saturated hydrocarbons with single C-C bonds only',
        emission_type='anthropogenic',
        subcategories={
            'c1_c4_light': VOCSubcategory(
                name='Light Alkanes (C1-C4)',
                description='Methane through butanes - highly volatile',
                compounds=['methane', 'ethane', 'propane', 'n-butane', 'isobutane'],
                typical_sources=['natural gas', 'vehicle exhaust', 'refining'],
                soa_potential='low'
            ),
            'c5_c10_medium': VOCSubcategory(
                name='Medium Alkanes (C5-C10)',
                description='Pentanes through decanes - gasoline range',
                compounds=['n_pentane', 'isopentane', 'neopentane', 'n_hexane',
                          '2_methylpentane', '3_methylpentane', '22_dimethylbutane',
                          '23_dimethylbutane', 'cyclohexane', 'methylcyclohexane',
                          'n_heptane', 'n_octane', 'isooctane', 'n_nonane', 'n_decane'],
                typical_sources=['gasoline evaporation', 'solvents', 'vehicle exhaust'],
                soa_potential='medium'
            ),
            'c11_c15_heavy': VOCSubcategory(
                name='Heavy Alkanes (C11-C15)',
                description='Undecanes through pentadecanes - diesel range',
                compounds=['n-undecane', 'n-dodecane', 'n-tridecane',
                          'n-tetradecane', 'n-pentadecane'],
                typical_sources=['diesel exhaust', 'cooking emissions', 'asphalt'],
                soa_potential='high'
            ),
            'cycloalkanes': VOCSubcategory(
                name='Cycloalkanes',
                description='Ring-containing saturated hydrocarbons',
                compounds=['cyclopentane', 'cyclohexane', 'methylcyclopentane',
                          'methylcyclohexane', 'ethylcyclohexane'],
                typical_sources=['gasoline', 'naphtha'],
                soa_potential='medium'
            ),
        }
    ),

    # =========================================================================
    # ALKENES (Unsaturated Hydrocarbons)
    # =========================================================================
    'alkenes': VOCCategory(
        name='Alkenes',
        description='Unsaturated hydrocarbons with C=C double bonds',
        emission_type='both',
        subcategories={
            'simple_alkenes': VOCSubcategory(
                name='Simple Alkenes',
                description='Non-branched alkenes',
                compounds=['ethene', 'propene', '1-butene', '2-butene',
                          'cis-2-butene', 'trans-2-butene', '1-pentene',
                          '1-hexene', '1-heptene', '1-octene'],
                typical_sources=['vehicle exhaust', 'petrochemical industry'],
                soa_potential='medium'
            ),
            'branched_alkenes': VOCSubcategory(
                name='Branched Alkenes',
                description='Methylated and branched alkenes',
                compounds=['isobutene', '2-methyl-1-butene', '2-methyl-2-butene',
                          '3-methyl-1-butene', '2-methylpropene'],
                typical_sources=['gasoline evaporation', 'chemical manufacturing'],
                soa_potential='medium'
            ),
            'dienes': VOCSubcategory(
                name='Dienes',
                description='Alkenes with two double bonds',
                compounds=['13_butadiene', '13_pentadiene', '14_pentadiene',
                          '23_dimethyl_13_butadiene', 'cyclopentadiene'],
                typical_sources=['vehicle exhaust', 'industrial processes'],
                soa_potential='high'
            ),
            'isoprene': VOCSubcategory(
                name='Isoprene',
                description='2-methyl-1,3-butadiene - most abundant biogenic VOC',
                compounds=['isoprene'],
                typical_sources=['deciduous trees', 'vegetation'],
                soa_potential='medium'
            ),
        }
    ),

    # =========================================================================
    # AROMATICS
    # =========================================================================
    'aromatics': VOCCategory(
        name='Aromatics',
        description='Benzene ring-containing compounds (BTEX and others)',
        emission_type='anthropogenic',
        subcategories={
            'btex': VOCSubcategory(
                name='BTEX',
                description='Benzene, Toluene, Ethylbenzene, Xylenes',
                compounds=['benzene', 'toluene', 'ethylbenzene',
                          'o-xylene', 'm-xylene', 'p-xylene'],
                typical_sources=['gasoline', 'solvents', 'vehicle exhaust'],
                soa_potential='high'
            ),
            'c9_aromatics': VOCSubcategory(
                name='C9 Aromatics',
                description='Trimethylbenzenes and related compounds',
                compounds=['123_trimethylbenzene', '124_trimethylbenzene',
                          '135_trimethylbenzene', 'propylbenzene',
                          'isopropylbenzene', 'indane', 'indene'],
                typical_sources=['gasoline', 'solvents', 'industrial'],
                soa_potential='very_high'
            ),
            'c10_aromatics': VOCSubcategory(
                name='C10+ Aromatics',
                description='Larger aromatics including naphthalene',
                compounds=['naphthalene', 'butylbenzene', 'diethylbenzene',
                          '1234_tetramethylbenzene', 'tetralin'],
                typical_sources=['diesel', 'coal tar', 'industrial'],
                soa_potential='very_high'
            ),
            'phenolics': VOCSubcategory(
                name='Phenolics',
                description='Hydroxylated aromatics',
                compounds=['phenol', 'cresols', 'catechol', 'guaiacol',
                          'syringol', 'methoxyphenols'],
                typical_sources=['biomass burning', 'wood combustion'],
                soa_potential='very_high'
            ),
            'styrenes': VOCSubcategory(
                name='Styrenes',
                description='Vinyl-substituted benzenes',
                compounds=['styrene', 'alpha-methylstyrene'],
                typical_sources=['polymer production', 'vehicle exhaust'],
                soa_potential='high'
            ),
        }
    ),

    # =========================================================================
    # MONOTERPENES (C10)
    # =========================================================================
    'monoterpenes': VOCCategory(
        name='Monoterpenes',
        description='C10 isoprenoids from vegetation (C10H16)',
        emission_type='biogenic',
        subcategories={
            'bicyclic': VOCSubcategory(
                name='Bicyclic Monoterpenes',
                description='Two-ring terpene structures',
                compounds=['alpha-pinene', 'beta-pinene', '3-carene',
                          'camphene', 'sabinene', 'alpha-thujene'],
                typical_sources=['coniferous trees', 'pine forests'],
                soa_potential='very_high'
            ),
            'monocyclic': VOCSubcategory(
                name='Monocyclic Monoterpenes',
                description='Single-ring terpene structures',
                compounds=['limonene', 'alpha-terpinene', 'gamma-terpinene',
                          'terpinolene', 'alpha-phellandrene', 'beta-phellandrene',
                          'p-cymene'],
                typical_sources=['citrus', 'cleaning products', 'broad-leaf trees'],
                soa_potential='very_high'
            ),
            'acyclic': VOCSubcategory(
                name='Acyclic Monoterpenes',
                description='Open-chain monoterpenes',
                compounds=['myrcene', 'beta-ocimene', 'alpha-ocimene',
                          'linalool', 'geraniol', 'citronellol'],
                typical_sources=['flowers', 'herbs', 'essential oils'],
                soa_potential='high'
            ),
        }
    ),

    # =========================================================================
    # SESQUITERPENES (C15)
    # =========================================================================
    'sesquiterpenes': VOCCategory(
        name='Sesquiterpenes',
        description='C15 isoprenoids from vegetation (C15H24)',
        emission_type='biogenic',
        subcategories={
            'major_sesquiterpenes': VOCSubcategory(
                name='Major Sesquiterpenes',
                description='Most abundant atmospheric sesquiterpenes',
                compounds=['beta-caryophyllene', 'alpha-humulene',
                          'beta-farnesene', 'alpha-farnesene',
                          'longifolene', 'alpha-cedrene'],
                typical_sources=['trees', 'agricultural crops', 'flowers'],
                soa_potential='very_high'
            ),
            'minor_sesquiterpenes': VOCSubcategory(
                name='Minor Sesquiterpenes',
                description='Less abundant sesquiterpenes',
                compounds=['alpha-bisabolene', 'beta-bisabolene',
                          'aromadendene', 'valencene', 'nerolidol'],
                typical_sources=['essential oils', 'specific plants'],
                soa_potential='very_high'
            ),
        }
    ),

    # =========================================================================
    # OXYGENATED VOCs (OVOCs)
    # =========================================================================
    'oxygenated': VOCCategory(
        name='Oxygenated VOCs',
        description='VOCs containing oxygen functional groups',
        emission_type='both',
        subcategories={
            'aldehydes': VOCSubcategory(
                name='Aldehydes',
                description='Compounds with -CHO group',
                compounds=['formaldehyde', 'acetaldehyde', 'propanal',
                          'butanal', 'pentanal', 'hexanal',
                          'acrolein', 'crotonaldehyde', 'methacrolein'],
                typical_sources=['vehicle exhaust', 'oxidation products', 'cooking'],
                soa_potential='medium'
            ),
            'ketones': VOCSubcategory(
                name='Ketones',
                description='Compounds with C=O group (not aldehyde)',
                compounds=['acetone', 'methyl_ethyl_ketone', 'methyl_vinyl_ketone',
                          '2-pentanone', '3-pentanone', 'methyl_isobutyl_ketone',
                          'cyclohexanone'],
                typical_sources=['solvents', 'oxidation products'],
                soa_potential='low'
            ),
            'alcohols': VOCSubcategory(
                name='Alcohols',
                description='Compounds with -OH group',
                compounds=['methanol', 'ethanol', 'isopropanol',
                          '1-butanol', '2-butanol', '1-propanol',
                          '2-methyl-1-propanol', '1-pentanol'],
                typical_sources=['solvents', 'biofuel', 'fermentation'],
                soa_potential='low'
            ),
            'carboxylic_acids': VOCSubcategory(
                name='Carboxylic Acids',
                description='Compounds with -COOH group',
                compounds=['formic_acid', 'acetic_acid', 'propionic_acid',
                          'butyric_acid', 'valeric_acid', 'hexanoic_acid'],
                typical_sources=['biomass burning', 'oxidation products'],
                soa_potential='high'
            ),
            'esters': VOCSubcategory(
                name='Esters',
                description='Compounds with -COO- linkage',
                compounds=['methyl_formate', 'methyl_acetate', 'ethyl_acetate',
                          'butyl_acetate', 'isoamyl_acetate'],
                typical_sources=['solvents', 'fragrances', 'fruits'],
                soa_potential='low'
            ),
            'ethers': VOCSubcategory(
                name='Ethers',
                description='Compounds with C-O-C linkage',
                compounds=['dimethyl_ether', 'diethyl_ether',
                          'methyl_tert_butyl_ether', 'ethyl_tert_butyl_ether',
                          'tetrahydrofuran', '2-methylfuran'],
                typical_sources=['fuel additives', 'solvents'],
                soa_potential='medium'
            ),
        }
    ),

    # =========================================================================
    # TERPENE OXIDATION PRODUCTS
    # =========================================================================
    'terpene_oxidation': VOCCategory(
        name='Terpene Oxidation Products',
        description='Products from atmospheric oxidation of terpenes',
        emission_type='biogenic',  # Secondary formation
        subcategories={
            'pinene_products': VOCSubcategory(
                name='Pinene Oxidation Products',
                description='Products from alpha- and beta-pinene oxidation',
                compounds=['pinonaldehyde', 'pinonic_acid', 'pinic_acid',
                          'norpinonaldehyde', 'nopinone', '10-hydroxypinonic_acid',
                          'MBTCA', 'terpenylic_acid'],
                typical_sources=['atmospheric oxidation of pinenes'],
                soa_potential='very_high'
            ),
            'limonene_products': VOCSubcategory(
                name='Limonene Oxidation Products',
                description='Products from limonene oxidation',
                compounds=['limonaldehyde', 'limonic_acid', 'limonalic_acid',
                          'ketolimonaldehyde', 'keto-limononaldehyde'],
                typical_sources=['atmospheric oxidation of limonene'],
                soa_potential='very_high'
            ),
            'isoprene_products': VOCSubcategory(
                name='Isoprene Oxidation Products',
                description='Products from isoprene oxidation',
                compounds=['methyl_vinyl_ketone', 'methacrolein',
                          'isoprene_epoxides', 'HMML', 'MAE',
                          'hydroxymethyl-methyl-alpha-lactone'],
                typical_sources=['atmospheric oxidation of isoprene'],
                soa_potential='medium'
            ),
            'sesquiterpene_products': VOCSubcategory(
                name='Sesquiterpene Oxidation Products',
                description='Products from sesquiterpene oxidation',
                compounds=['caryophyllonic_acid', 'caryophyllinic_acid',
                          'beta-nocaryophyllone'],
                typical_sources=['atmospheric oxidation of sesquiterpenes'],
                soa_potential='very_high'
            ),
        }
    ),
}


# ==============================================================================
# Helper Functions
# ==============================================================================

def get_all_compounds() -> List[str]:
    """Get list of all compound names across all categories."""
    compounds = []
    for category in VOC_CATEGORIES.values():
        for subcategory in category.subcategories.values():
            compounds.extend(subcategory.compounds)
    return sorted(set(compounds))


def get_category_compounds(category_name: str) -> List[str]:
    """Get all compounds in a category."""
    if category_name not in VOC_CATEGORIES:
        return []

    category = VOC_CATEGORIES[category_name]
    compounds = []
    for subcategory in category.subcategories.values():
        compounds.extend(subcategory.compounds)
    return sorted(set(compounds))


def get_subcategory_compounds(category_name: str, subcategory_name: str) -> List[str]:
    """Get all compounds in a specific subcategory."""
    if category_name not in VOC_CATEGORIES:
        return []

    category = VOC_CATEGORIES[category_name]
    if subcategory_name not in category.subcategories:
        return []

    return category.subcategories[subcategory_name].compounds.copy()


def get_compound_category(compound_name: str) -> Optional[Dict[str, str]]:
    """Find the category and subcategory for a compound."""
    compound_lower = compound_name.lower().replace('-', '').replace('_', '').replace(' ', '')

    for cat_name, category in VOC_CATEGORIES.items():
        for subcat_name, subcategory in category.subcategories.items():
            for comp in subcategory.compounds:
                comp_lower = comp.lower().replace('-', '').replace('_', '').replace(' ', '')
                if comp_lower == compound_lower:
                    return {
                        'category': cat_name,
                        'subcategory': subcat_name,
                        'emission_type': category.emission_type,
                        'soa_potential': subcategory.soa_potential
                    }
    return None


def get_high_soa_compounds() -> List[str]:
    """Get compounds with high or very high SOA potential."""
    compounds = []
    for category in VOC_CATEGORIES.values():
        for subcategory in category.subcategories.values():
            if subcategory.soa_potential in ('high', 'very_high'):
                compounds.extend(subcategory.compounds)
    return sorted(set(compounds))


def get_biogenic_compounds() -> List[str]:
    """Get all biogenic VOCs."""
    compounds = []
    for category in VOC_CATEGORIES.values():
        if category.emission_type in ('biogenic', 'both'):
            for subcategory in category.subcategories.values():
                compounds.extend(subcategory.compounds)
    return sorted(set(compounds))


def get_anthropogenic_compounds() -> List[str]:
    """Get all anthropogenic VOCs."""
    compounds = []
    for category in VOC_CATEGORIES.values():
        if category.emission_type in ('anthropogenic', 'both'):
            for subcategory in category.subcategories.values():
                compounds.extend(subcategory.compounds)
    return sorted(set(compounds))


def get_category_statistics() -> Dict[str, Dict[str, int]]:
    """Get statistics about VOC categories."""
    stats = {}
    for cat_name, category in VOC_CATEGORIES.items():
        total_compounds = 0
        subcats = len(category.subcategories)
        for subcategory in category.subcategories.values():
            total_compounds += len(subcategory.compounds)

        stats[cat_name] = {
            'subcategories': subcats,
            'compounds': total_compounds,
            'emission_type': category.emission_type
        }
    return stats


# ==============================================================================
# Dropdown Options Generator
# ==============================================================================

def generate_dropdown_options(group_by: str = 'category') -> List[Dict[str, str]]:
    """
    Generate dropdown options for UI, optionally grouped.

    Args:
        group_by: 'category', 'subcategory', 'soa_potential', or 'emission_type'

    Returns:
        List of dicts with 'label' and 'value' keys
    """
    try:
        # Import compound database
        from .compound_database import CompoundDatabase
        db = CompoundDatabase()

        if group_by == 'category':
            options = []
            for cat_name, category in VOC_CATEGORIES.items():
                # Add category header
                options.append({
                    'label': f"── {category.name} ──",
                    'value': '',
                    'disabled': True
                })
                # Add compounds in category
                for subcat in category.subcategories.values():
                    for comp_name in sorted(subcat.compounds):
                        compound = db.get_by_name(comp_name)
                        if compound:
                            options.append({
                                'label': compound.name,
                                'value': compound.name.lower().replace(' ', '-')
                            })
            return options

        elif group_by == 'emission_type':
            options = []
            for emission_type in ['biogenic', 'anthropogenic', 'both']:
                compounds = []
                for category in VOC_CATEGORIES.values():
                    if category.emission_type == emission_type:
                        for subcat in category.subcategories.values():
                            compounds.extend(subcat.compounds)

                if compounds:
                    options.append({
                        'label': f"── {emission_type.title()} ──",
                        'value': '',
                        'disabled': True
                    })
                    for comp_name in sorted(set(compounds)):
                        compound = db.get_by_name(comp_name)
                        if compound:
                            options.append({
                                'label': compound.name,
                                'value': compound.name.lower().replace(' ', '-')
                            })
            return options

        else:
            # Simple alphabetical list
            return [
                {
                    'label': compound.name,
                    'value': compound.name.lower().replace(' ', '-')
                }
                for compound in sorted(db.compounds.values(), key=lambda c: c.name)
            ]

    except ImportError:
        # Fallback if compound database not available
        options = []
        for cat_name, category in VOC_CATEGORIES.items():
            for subcat in category.subcategories.values():
                for comp_name in sorted(subcat.compounds):
                    options.append({
                        'label': comp_name,
                        'value': comp_name.lower().replace(' ', '-')
                    })
        return options


# ==============================================================================
# Module Test
# ==============================================================================

if __name__ == "__main__":
    print("VOC Categories Summary")
    print("=" * 60)

    stats = get_category_statistics()
    total_compounds = sum(s['compounds'] for s in stats.values())

    for cat_name, stat in stats.items():
        print(f"\n{cat_name.upper()}:")
        print(f"  Subcategories: {stat['subcategories']}")
        print(f"  Compounds: {stat['compounds']}")
        print(f"  Emission type: {stat['emission_type']}")

    print(f"\n{'=' * 60}")
    print(f"Total unique compounds: {len(get_all_compounds())}")
    print(f"Biogenic VOCs: {len(get_biogenic_compounds())}")
    print(f"Anthropogenic VOCs: {len(get_anthropogenic_compounds())}")
    print(f"High SOA potential: {len(get_high_soa_compounds())}")
