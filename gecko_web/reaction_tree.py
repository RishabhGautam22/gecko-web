"""
GECKO-A Reaction Tree Parser and SMILES Converter

This module handles parsing of GECKO-A mechanism files and conversion of
GECKO formula notation to standard SMILES for visualization.

Author: Deeksha Sharma
"""

import os
import re
import json
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ==============================================================================
# Comprehensive KNOWN_SPECIES Database
# ==============================================================================
# This database is the AUTHORITATIVE source for SMILES conversions.
# It maps GECKO codes, formulas, and common names to verified SMILES strings.
# Entries have been validated with RDKit for correctness.

KNOWN_SPECIES = {
    # =========================================================================
    # INORGANICS & SMALL MOLECULES
    # =========================================================================
    'H2O': 'O', 'WATER': 'O',
    'H2O2': 'OO', 'HYDROGENPEROXIDE': 'OO',
    'HNO3': 'O[N+](=O)[O-]', 'NITRICACID': 'O[N+](=O)[O-]',
    'HNO2': 'ON=O', 'NITROUSACID': 'ON=O',
    'N2O5': 'O=[N+]([O-])O[N+](=O)[O-]',
    'CH2O': 'C=O', 'HCHO': 'C=O', 'FORMALDEHYDE': 'C=O',
    'CO': '[C-]#[O+]', 'CARBONMONOXIDE': '[C-]#[O+]',
    'CO2': 'O=C=O', 'CARBONDIOXIDE': 'O=C=O',
    'NO': '[N]=O', 'NITRICOXIDE': '[N]=O',
    'NO2': '[N+](=O)[O-]', 'NITROGENDIOXIDE': '[N+](=O)[O-]',
    'NO3': '[O][N+](=O)[O-]', 'NITRATE': '[O][N+](=O)[O-]',
    'O3': '[O-][O+]=O', 'OZONE': '[O-][O+]=O',
    'HO': '[OH]', 'OH': '[OH]', 'HYDROXYL': '[OH]',
    'HO2': 'O[O]', 'HYDROPEROXYL': 'O[O]',
    'CH4': 'C', 'METHANE': 'C',
    'HCOOH': 'OC=O', 'FORMICACID': 'OC=O', 'FORMIC': 'OC=O',
    'CH3OH': 'CO', 'METHANOL': 'CO',
    'SO2': 'O=S=O', 'SULFURDIOXIDE': 'O=S=O',
    'H2SO4': 'OS(=O)(=O)O', 'SULFURICACID': 'OS(=O)(=O)O',
    'H2': '[H][H]', 'HYDROGEN': '[H][H]',
    'O2': 'O=O', 'OXYGEN': 'O=O',
    'N2': 'N#N', 'NITROGEN': 'N#N',

    # =========================================================================
    # C1 COMPOUNDS
    # =========================================================================
    'CH3OOH': 'COO', 'METHYLHYDROPEROXIDE': 'COO', 'MHP': 'COO',
    'CH3ONO2': 'CO[N+](=O)[O-]', 'METHYLNITRATE': 'CO[N+](=O)[O-]',
    'CH3OO': 'CO[O]', 'CH3O2': 'CO[O]', 'METHYLPEROXY': 'CO[O]',
    'CH3O': 'C[O]', 'METHOXY': 'C[O]',
    'CH3NO3': 'CO[N+](=O)[O-]',

    # =========================================================================
    # C2 COMPOUNDS
    # =========================================================================
    'C2H6': 'CC', 'ETHANE': 'CC', 'C02000': 'CC',
    'C2H4': 'C=C', 'ETHENE': 'C=C', 'ETHYLENE': 'C=C',
    'C2H2': 'C#C', 'ETHYNE': 'C#C', 'ACETYLENE': 'C#C',
    'CH3CHO': 'CC=O', 'ACETALD': 'CC=O', 'ACETALDEHYDE': 'CC=O', 'ACD': 'CC=O',
    'GLYOX': 'O=CC=O', 'GLYOXAL': 'O=CC=O', 'GLY': 'O=CC=O',
    'C2H5OH': 'CCO', 'ETHANOL': 'CCO', 'ETH': 'CCO',
    'C2H5OOH': 'CCOO', 'ETHYLHYDROPEROXIDE': 'CCOO',
    'C2H5ONO2': 'CCO[N+](=O)[O-]', 'ETHYLNITRATE': 'CCO[N+](=O)[O-]',
    'C2H5OO': 'CCO[O]', 'ETHYLPEROXY': 'CCO[O]',
    'C2H5O': 'CC[O]', 'ETHOXY': 'CC[O]',
    'CH3COOH': 'CC(=O)O', 'ACETICACID': 'CC(=O)O', 'ACETIC': 'CC(=O)O',
    'PAN': 'CC(=O)OO[N+](=O)[O-]', 'PEROXYACETYLNITRATE': 'CC(=O)OO[N+](=O)[O-]',
    'PERACETIC': 'CC(=O)OO', 'PAA': 'CC(=O)OO',
    'CH3CO': 'CC(=O)[O]', 'ACETYL': 'CC(=O)[O]',
    'CH3CO3': 'CC(=O)O[O]', 'ACETYLPEROXY': 'CC(=O)O[O]',
    'OXALICACID': 'OC(=O)C(=O)O', 'OXALIC': 'OC(=O)C(=O)O',

    # =========================================================================
    # C3 COMPOUNDS
    # =========================================================================
    'C3H8': 'CCC', 'PROPANE': 'CCC', 'C03000': 'CCC',
    'C3H6': 'CC=C', 'PROPENE': 'CC=C', 'PROPYLENE': 'CC=C',
    'ACETONE': 'CC(=O)C', 'CH3COCH3': 'CC(=O)C', 'ACE': 'CC(=O)C',
    'PROPANAL': 'CCC=O', 'PROPIONALDEHYDE': 'CCC=O',
    'MGLYOX': 'CC(=O)C=O', 'METHYLGLYOXAL': 'CC(=O)C=O', 'MGLY': 'CC(=O)C=O',
    'HYAC': 'CC(=O)CO', 'HYDROXYACETONE': 'CC(=O)CO', 'HACET': 'CC(=O)CO',
    'ACROLEIN': 'C=CC=O', 'ACR': 'C=CC=O',
    'PROPANOL': 'CCCO', '1PROPANOL': 'CCCO',
    '2PROPANOL': 'CC(O)C', 'ISOPROPANOL': 'CC(O)C', 'IPA': 'CC(O)C',
    'PROPYLNITRATE': 'CCCO[N+](=O)[O-]',
    'PROPYLPEROXY': 'CCCO[O]',
    'LACTIC': 'CC(O)C(=O)O', 'LACTICACID': 'CC(O)C(=O)O',
    'PYRUVIC': 'CC(=O)C(=O)O', 'PYRUVICACID': 'CC(=O)C(=O)O',
    'PPAN': 'CCC(=O)OO[N+](=O)[O-]', 'PEROXYPROPIONYLNITRATE': 'CCC(=O)OO[N+](=O)[O-]',
    'PROPIONICACID': 'CCC(=O)O', 'PROPIONIC': 'CCC(=O)O',

    # =========================================================================
    # C4 COMPOUNDS
    # =========================================================================
    'C4H10': 'CCCC', 'BUTANE': 'CCCC', 'NBUTANE': 'CCCC', 'C04000': 'CCCC',
    'IBUTANE': 'CC(C)C', 'ISOBUTANE': 'CC(C)C', '2METHYLPROPANE': 'CC(C)C',
    '1BUTENE': 'CCC=C', 'BUTENE': 'CCC=C',
    '2BUTENE': 'CC=CC', 'CIS2BUTENE': 'C/C=C\\C', 'TRANS2BUTENE': 'C/C=C/C',
    'IBUTENE': 'CC(=C)C', 'ISOBUTENE': 'CC(=C)C', '2METHYLPROPENE': 'CC(=C)C',
    '13BUTADIENE': 'C=CC=C', 'BUTADIENE': 'C=CC=C', '1,3BUTADIENE': 'C=CC=C',
    'MVK': 'CC(=O)C=C', 'METHYLVINYLKETONE': 'CC(=O)C=C',
    'MACR': 'CC(=C)C=O', 'METHACROLEIN': 'CC(=C)C=O',
    'MEK': 'CCC(=O)C', 'METHYLETHYLKETONE': 'CCC(=O)C', 'BUTANONE': 'CCC(=O)C',
    'BUTANAL': 'CCCC=O', 'BUTYRALDEHYDE': 'CCCC=O',
    'ISOBUTANAL': 'CC(C)C=O', 'ISOBUTYRALDEHYDE': 'CC(C)C=O',
    'FURAN': 'c1ccoc1', 'FURANE': 'c1ccoc1',
    'THF': 'C1CCCO1', 'TETRAHYDROFURAN': 'C1CCCO1',
    '2METHYLFURAN': 'Cc1ccoc1', 'MFURAN': 'Cc1ccoc1',
    '3METHYLFURAN': 'Cc1ccoc1',
    'BUTYRICACID': 'CCCC(=O)O', 'BUTYRIC': 'CCCC(=O)O',
    'ISOBUTYRICACID': 'CC(C)C(=O)O', 'ISOBUTYRIC': 'CC(C)C(=O)O',
    'BUTANOL': 'CCCCO', '1BUTANOL': 'CCCCO',
    '2BUTANOL': 'CCC(O)C',
    'ISOBUTANOL': 'CC(C)CO',
    'TERTBUTANOL': 'CC(C)(C)O', 'TBUTANOL': 'CC(C)(C)O',
    'CROTONALDEHYDE': 'CC=CC=O',
    'MALEICACID': 'OC(=O)C=CC(=O)O', 'MALEIC': 'OC(=O)C=CC(=O)O',
    'FUMARICACID': 'OC(=O)C=CC(=O)O', 'FUMARIC': 'OC(=O)C=CC(=O)O',
    'SUCCINICACID': 'OC(=O)CCC(=O)O', 'SUCCINIC': 'OC(=O)CCC(=O)O',
    'MALIC': 'OC(CC(=O)O)C(=O)O', 'MALICACID': 'OC(CC(=O)O)C(=O)O',
    'TARTARIC': 'OC(C(O)C(=O)O)C(=O)O', 'TARTARICACID': 'OC(C(O)C(=O)O)C(=O)O',

    # =========================================================================
    # C5 COMPOUNDS - ISOPRENE AND DERIVATIVES
    # =========================================================================
    'C5H8': 'CC(=C)C=C', 'ISOP': 'CC(=C)C=C',
    'ISOPRE': 'CC(=C)C=C', 'ISOPRENE': 'CC(=C)C=C',
    'C5H12': 'CCCCC', 'PENTANE': 'CCCCC', 'NPENTANE': 'CCCCC', 'C05000': 'CCCCC',
    'IPENTANE': 'CC(C)CC', 'ISOPENTANE': 'CC(C)CC', '2METHYLBUTANE': 'CC(C)CC', 'IPENTA': 'CC(C)CC',
    # GECKO formula notation for isopentane: CH3CH(CH3)CH2CH3
    'CH3CH(CH3)CH2CH3': 'CC(C)CC',
    'NEOPENTANE': 'CC(C)(C)C', '2,2DIMETHYLPROPANE': 'CC(C)(C)C', 'NEOPEN': 'CC(C)(C)C',
    # GECKO formula notation for neopentane: CH3C(CH3)(CH3)CH3
    'CH3C(CH3)(CH3)CH3': 'CC(C)(C)C',
    '1PENTENE': 'CCCC=C',
    '2PENTENE': 'CCC=CC',
    '2METHYL1BUTENE': 'CCC(=C)C',
    '2METHYL2BUTENE': 'CC=C(C)C',
    '3METHYL1BUTENE': 'CC(C)C=C',
    'CYCLOPENTANE': 'C1CCCC1',
    'CYCLOPENTENE': 'C1CCC=C1',
    'CYCLOPENTADIENE': 'C1C=CC=C1', 'CPD': 'C1C=CC=C1',
    'METHYLCYCLOPENTANE': 'CC1CCCC1',
    # Isoprene oxidation products
    'ISOPAO2': 'CC(=C)C(O[O])C=C', 'ISOPA_OO': 'CC(=C)C(O[O])C=C',
    'ISOPBO2': 'CC(O[O])(C=C)C=C', 'ISOPB_OO': 'CC(O[O])(C=C)C=C',
    'ISOPCO2': 'CC(=C)C(C=C)O[O]', 'ISOPC_OO': 'CC(=C)C(C=C)O[O]',
    'ISOPDO2': 'CC(O[O])=CC=C', 'ISOPD_OO': 'CC(O[O])=CC=C',
    'ISOPAOH': 'CC(=C)C(O)C=C', 'ISOPA_OH': 'CC(=C)C(O)C=C',
    'ISOPBOH': 'CC(O)(C=C)C=C', 'ISOPB_OH': 'CC(O)(C=C)C=C',
    'ISOPCOH': 'CC(=C)C(C=C)O', 'ISOPC_OH': 'CC(=C)C(C=C)O',
    'ISOPDOH': 'CC(O)=CC=C', 'ISOPD_OH': 'CC(O)=CC=C',
    'ISOPAOOH': 'CC(=C)C(OO)C=C', 'ISOPA_OOH': 'CC(=C)C(OO)C=C',
    'ISOPBOOH': 'CC(OO)(C=C)C=C', 'ISOPB_OOH': 'CC(OO)(C=C)C=C',
    'ISOPCOOH': 'CC(=C)C(C=C)OO', 'ISOPC_OOH': 'CC(=C)C(C=C)OO',
    'ISOPDOOH': 'CC(OO)=CC=C', 'ISOPD_OOH': 'CC(OO)=CC=C',
    'ISOPANO3': 'CC(=C)C(O[N+](=O)[O-])C=C',
    'ISOPBNO3': 'CC(O[N+](=O)[O-])(C=C)C=C',
    'ISOPCNO3': 'CC(=C)C(C=C)O[N+](=O)[O-]',
    'ISOPDNO3': 'CC(O[N+](=O)[O-])=CC=C',
    'HMACR': 'CC(=O)C(C)=CO', 'HMETHACROLEIN': 'CC(=O)C(C)=CO',
    'HMPAN': 'OCC(=C)C(=O)OO[N+](=O)[O-]',
    'IEPOX': 'CC1(O)OC1C=C', 'ISOPRENEEPOXIDE': 'CC1(O)OC1C=C',
    'ISOPEPOX': 'CC1OC1C(O)C', 'ISOP_EPOX': 'CC1OC1C(O)C',
    'GLUTARIC': 'OC(=O)CCCC(=O)O', 'GLUTARICACID': 'OC(=O)CCCC(=O)O',
    'VALERICACID': 'CCCCC(=O)O', 'VALERIC': 'CCCCC(=O)O',
    'ISOVALERIC': 'CC(C)CC(=O)O', 'ISOVALERICACID': 'CC(C)CC(=O)O',
    'PENTANAL': 'CCCCC=O', 'VALERALDEHYDE': 'CCCCC=O',
    '2PENTANONE': 'CCCC(=O)C',
    '3PENTANONE': 'CCC(=O)CC',
    'PENTANOL': 'CCCCCO', '1PENTANOL': 'CCCCCO',
    '2PENTANOL': 'CCCC(O)C',
    '3PENTANOL': 'CCC(O)CC',

    # =========================================================================
    # C6 COMPOUNDS
    # =========================================================================
    'C6H14': 'CCCCCC', 'HEXANE': 'CCCCCC', 'NHEXANE': 'CCCCCC', 'C06000': 'CCCCCC',
    '2METHYLPENTANE': 'CCCC(C)C', 'ISOHEXANE': 'CCCC(C)C',
    '3METHYLPENTANE': 'CCC(C)CC',
    '2,2DIMETHYLBUTANE': 'CCC(C)(C)C', 'NEOHEXANE': 'CCC(C)(C)C',
    '2,3DIMETHYLBUTANE': 'CC(C)C(C)C',
    'CYCLOHEXANE': 'C1CCCCC1', 'CHEX': 'C1CCCCC1',
    'CYCLOHEXENE': 'C1CCC=CC1', 'CHEXENE': 'C1CCC=CC1',
    'METHYLCYCLOHEXANE': 'CC1CCCCC1', 'MCHEX': 'CC1CCCCC1',
    '1,3CYCLOHEXADIENE': 'C1C=CC=CC1',
    '1,4CYCLOHEXADIENE': 'C1C=CCC=C1',
    'BENZENE': 'c1ccccc1', 'BENZ': 'c1ccccc1', 'C6H6': 'c1ccccc1',
    'R06000': 'c1ccccc1',  # GECKO internal code for benzene (C6 aromatic)
    '1HEXENE': 'CCCCC=C',
    '2HEXENE': 'CCCC=CC',
    '3HEXENE': 'CCC=CCC',
    '2,4HEXADIENE': 'CC=CC=CC',
    '1,3HEXADIENE': 'C=CC=CCC',
    '1,5HEXADIENE': 'C=CCCC=C',
    'HEXANAL': 'CCCCCC=O', 'CAPROALDEHYDE': 'CCCCCC=O',
    '2HEXANONE': 'CCCCC(=O)C',
    '3HEXANONE': 'CCCC(=O)CC',
    'ADIPICACID': 'OC(=O)CCCCC(=O)O', 'ADIPIC': 'OC(=O)CCCCC(=O)O',
    'CAPROICACID': 'CCCCCC(=O)O', 'HEXANOICACID': 'CCCCCC(=O)O',
    'HEXANOL': 'CCCCCCO', '1HEXANOL': 'CCCCCCO',
    '2HEXANOL': 'CCCCC(O)C',
    'CYCLOHEXANOL': 'OC1CCCCC1',
    'CYCLOHEXANONE': 'O=C1CCCCC1',
    'SORBITOL': 'OCC(O)C(O)C(O)C(O)CO',
    'PHENOL': 'Oc1ccccc1', 'PHENYL': 'Oc1ccccc1',
    'CATECHOL': 'Oc1ccccc1O', '1,2BENZENEDIOL': 'Oc1ccccc1O',
    'RESORCINOL': 'Oc1cccc(O)c1', '1,3BENZENEDIOL': 'Oc1cccc(O)c1',
    'HYDROQUINONE': 'Oc1ccc(O)cc1', '1,4BENZENEDIOL': 'Oc1ccc(O)cc1',
    'PYROGALLOL': 'Oc1cccc(O)c1O', '1,2,3BENZENETRIOL': 'Oc1cccc(O)c1O',

    # =========================================================================
    # C7 COMPOUNDS - TOLUENE AND DERIVATIVES
    # =========================================================================
    'C7H16': 'CCCCCCC', 'HEPTANE': 'CCCCCCC', 'NHEPTANE': 'CCCCCCC', 'C07000': 'CCCCCCC',
    'TOLUENE': 'Cc1ccccc1', 'TOLUEN': 'Cc1ccccc1', 'TOL': 'Cc1ccccc1',
    'METHYLBENZENE': 'Cc1ccccc1', 'C6H5CH3': 'Cc1ccccc1', 'C7H8': 'Cc1ccccc1',
    'R07000': 'Cc1ccccc1',  # GECKO internal code for toluene
    # Toluene oxidation products (GECKO internal codes)
    '2R7000': 'Cc1ccccc1O[O]',  # Benzyl peroxy radical (corrected)
    '1R7000': 'Cc1ccccc1[O]',   # Benzyloxy radical (corrected)
    'DR7000': 'O=Cc1ccccc1',    # Benzaldehyde
    'NR7000': 'Cc1ccccc1O[N+](=O)[O-]',  # Benzyl nitrate (corrected)
    'HR7000': 'Cc1ccccc1OO',    # Benzyl hydroperoxide (corrected)
    'OR7000': 'Cc1ccccc1O',     # Cresol (ortho)
    'OR7001': 'OCc1ccccc1',     # Benzyl alcohol
    '2R7001': 'Cc1ccc(O[O])cc1', # Cresyl peroxy
    '1R7001': 'Cc1ccc([O])cc1', # Cresoxy radical
    'TR7000': '[O]Cc1ccccc1',   # Benzyl radical
    'PR7000': 'OOCc1ccccc1',    # Benzyl hydroperoxide alt
    'BENZALDEHYDE': 'O=Cc1ccccc1', 'BZALDEHYDE': 'O=Cc1ccccc1', 'BENZALD': 'O=Cc1ccccc1',
    'CRESOL': 'Cc1ccccc1O', 'OCRESOL': 'Cc1ccccc1O', 'ORTHOCRESOL': 'Cc1ccccc1O',
    'MCRESOL': 'Cc1cccc(O)c1', 'METACRESOL': 'Cc1cccc(O)c1',
    'PCRESOL': 'Cc1ccc(O)cc1', 'PARACRESOL': 'Cc1ccc(O)cc1',
    'BENZYLALC': 'OCc1ccccc1', 'BENZYLALCOHOL': 'OCc1ccccc1', 'BZALC': 'OCc1ccccc1',
    'BENZOICACID': 'OC(=O)c1ccccc1', 'BENZOIC': 'OC(=O)c1ccccc1',
    'BENZYLHYDROPEROXIDE': 'OOCc1ccccc1',
    # More toluene oxidation products
    'NO7000': 'Cc1ccc(O[N+](=O)[O-])c(O)c1',  # Nitrocresol
    'VO7000': 'Cc1ccc([N+](=O)[O-])c(O)c1',  # Nitrocresol isomer
    '2T700O': 'OC(C=CC=CC(CO[N+](=O)[O-])O[O])=O',  # Ring-opened product
    # GECKO codes for cresols
    'R70100': 'Cc1ccccc1O',     # o-Cresol (2-methylphenol)
    'R70200': 'Cc1cccc(O)c1',   # m-Cresol (3-methylphenol)
    'R70300': 'Cc1ccc(O)cc1',   # p-Cresol (4-methylphenol)
    # Ring-opening products
    'MGLYO7': 'CC(=O)C=CC=O',   # Methylglyoxal-like ring-open
    'GLYO7': 'O=CC=CC=O',       # Glyoxal-like ring-open
    'HEPTANAL': 'CCCCCCC=O',
    '2HEPTANONE': 'CCCCCC(=O)C',
    '3HEPTANONE': 'CCCCC(=O)CC',
    '4HEPTANONE': 'CCCC(=O)CCC',
    'HEPTANOL': 'CCCCCCCO', '1HEPTANOL': 'CCCCCCCO',
    'HEPTANOICACID': 'CCCCCCC(=O)O',
    'METHYLCYCLOHEXENE': 'CC1CCC=CC1',
    'PIMELICACID': 'OC(=O)CCCCCC(=O)O', 'PIMELIC': 'OC(=O)CCCCCC(=O)O',
    'SALICYLALDEHYDE': 'O=Cc1ccccc1O',
    'SALICYLICACID': 'OC(=O)c1ccccc1O', 'SALICYLIC': 'OC(=O)c1ccccc1O',
    'ANISOLE': 'COc1ccccc1', 'METHOXYBENZENE': 'COc1ccccc1',
    'GUAIACOL': 'COc1ccccc1O', '2METHOXYPHENOL': 'COc1ccccc1O',

    # =========================================================================
    # C8 COMPOUNDS - XYLENES AND DERIVATIVES
    # =========================================================================
    'C8H18': 'CCCCCCCC', 'OCTANE': 'CCCCCCCC', 'NOCTANE': 'CCCCCCCC', 'C08000': 'CCCCCCCC',
    'ISOOCTANE': 'CC(C)CC(C)(C)C', '2,2,4TRIMETHYLPENTANE': 'CC(C)CC(C)(C)C',
    'OXYLENE': 'Cc1ccccc1C', 'OXYL': 'Cc1ccccc1C', '1,2DIMETHYLBENZENE': 'Cc1ccccc1C',
    'R08100': 'Cc1ccccc1C',  # GECKO code for o-xylene
    'MXYLENE': 'Cc1cccc(C)c1', 'MXYL': 'Cc1cccc(C)c1', '1,3DIMETHYLBENZENE': 'Cc1cccc(C)c1',
    'R08200': 'Cc1cccc(C)c1',  # GECKO code for m-xylene
    'PXYLENE': 'Cc1ccc(C)cc1', 'PXYL': 'Cc1ccc(C)cc1', '1,4DIMETHYLBENZENE': 'Cc1ccc(C)cc1',
    'R08300': 'Cc1ccc(C)cc1',  # GECKO code for p-xylene
    'XYLENE': 'Cc1ccc(C)cc1', 'XYL': 'Cc1ccc(C)cc1',
    'ETHYLBENZENE': 'CCc1ccccc1', 'EBENZ': 'CCc1ccccc1', 'ETHBENZ': 'CCc1ccccc1',
    'R08000': 'CCc1ccccc1',  # GECKO internal code for ethylbenzene (C8 aromatic)
    'STYRENE': 'C=Cc1ccccc1', 'VINYLBENZENE': 'C=Cc1ccccc1', 'PHENYLETHENE': 'C=Cc1ccccc1',
    # Xylene oxidation products
    'DIMETHYLPHENOL': 'Cc1cc(O)ccc1C', '2,4DIMETHYLPHENOL': 'Cc1cc(O)ccc1C',
    'DIMETHYLBENZALDEHYDE': 'Cc1ccc(C=O)cc1C',
    # Common C8 organics
    'OCTANAL': 'CCCCCCCC=O', 'CAPRYLALDEHYDE': 'CCCCCCCC=O',
    '2OCTANONE': 'CCCCCCC(=O)C',
    '3OCTANONE': 'CCCCCC(=O)CC',
    'OCTANOL': 'CCCCCCCCO', '1OCTANOL': 'CCCCCCCCO',
    '2OCTANOL': 'CCCCCCC(O)C',
    'OCTANOICACID': 'CCCCCCCC(=O)O', 'CAPRYLICACID': 'CCCCCCCC(=O)O',
    'SUBERICACID': 'OC(=O)CCCCCCC(=O)O', 'SUBERIC': 'OC(=O)CCCCCCC(=O)O',
    '2ETHYLTOLUENE': 'CCc1ccccc1C',
    '3ETHYLTOLUENE': 'CCc1cccc(C)c1',
    '4ETHYLTOLUENE': 'CCc1ccc(C)cc1',
    'CUMENE': 'CC(C)c1ccccc1', 'ISOPROPYLBENZENE': 'CC(C)c1ccccc1',
    'PROPYLBENZENE': 'CCCc1ccccc1', 'NPROPYLBENZENE': 'CCCc1ccccc1',
    '1METHYLSTYRENE': 'CC(=C)c1ccccc1', 'ALPHAMETHYLSTYRENE': 'CC(=C)c1ccccc1',

    # =========================================================================
    # C9 COMPOUNDS
    # =========================================================================
    'C9H20': 'CCCCCCCCC', 'NONANE': 'CCCCCCCCC', 'NNONANE': 'CCCCCCCCC', 'C09000': 'CCCCCCCCC',
    'TRIMETHYLBENZENE': 'Cc1cc(C)cc(C)c1', '1,3,5TRIMETHYLBENZENE': 'Cc1cc(C)cc(C)c1',
    'R09000': 'Cc1cc(C)cc(C)c1',  # GECKO code for 1,3,5-trimethylbenzene (mesitylene)
    '1,2,3TRIMETHYLBENZENE': 'Cc1cccc(C)c1C', 'HEMIMELLITENE': 'Cc1cccc(C)c1C',
    'R09100': 'Cc1cccc(C)c1C',  # GECKO code for 1,2,3-trimethylbenzene
    '1,2,4TRIMETHYLBENZENE': 'Cc1ccc(C)c(C)c1', 'PSEUDOCUMENE': 'Cc1ccc(C)c(C)c1',
    'R09200': 'Cc1ccc(C)c(C)c1',  # GECKO code for 1,2,4-trimethylbenzene
    'MESITYLENE': 'Cc1cc(C)cc(C)c1',
    'INDANE': 'C1Cc2ccccc2C1',
    'INDENE': 'C1=Cc2ccccc2C1',
    'NONANAL': 'CCCCCCCCC=O',
    'NONANOL': 'CCCCCCCCCO', '1NONANOL': 'CCCCCCCCCO',
    'NONANOICACID': 'CCCCCCCCC(=O)O', 'PELARGONIC': 'CCCCCCCCC(=O)O',
    'AZELAIC': 'OC(=O)CCCCCCCC(=O)O', 'AZELAICACID': 'OC(=O)CCCCCCCC(=O)O',

    # =========================================================================
    # C10 COMPOUNDS - TERPENES AND DERIVATIVES
    # =========================================================================
    'C10H22': 'CCCCCCCCCC', 'DECANE': 'CCCCCCCCCC', 'NDECANE': 'CCCCCCCCCC', 'C10000': 'CCCCCCCCCC',
    # Monoterpenes (C10H16)
    'APINENE': 'CC1=CCC2CC1C2(C)C', 'APINEN': 'CC1=CCC2CC1C2(C)C',
    'ALPHAPINENE': 'CC1=CCC2CC1C2(C)C', 'A-PINENE': 'CC1=CCC2CC1C2(C)C',
    'ALPHA-PINENE': 'CC1=CCC2CC1C2(C)C', 'ALPINENE': 'CC1=CCC2CC1C2(C)C',
    'BPINENE': 'CC1(C)C2CCC(=C)C1C2', 'BPINEN': 'CC1(C)C2CCC(=C)C1C2',
    'BETAPINENE': 'CC1(C)C2CCC(=C)C1C2', 'B-PINENE': 'CC1(C)C2CCC(=C)C1C2',
    'BETA-PINENE': 'CC1(C)C2CCC(=C)C1C2',
    'LIMONENE': 'CC1=CCC(CC1)C(=C)C', 'LIMONE': 'CC1=CCC(CC1)C(=C)C',
    'D-LIMONENE': 'CC1=CCC(CC1)C(=C)C', 'DLIMONENE': 'CC1=CCC(CC1)C(=C)C',
    'MYRCENE': 'CC(=C)CCCC(=C)C=C', 'MYRCEN': 'CC(=C)CCCC(=C)C=C',
    'BETAMYRCENE': 'CC(=C)CCCC(=C)C=C', 'B-MYRCENE': 'CC(=C)CCCC(=C)C=C',
    'OCIMENE': 'CC(=C)C=CCC(=C)C=C', 'OCIMEN': 'CC(=C)C=CCC(=C)C=C',
    'BETAOCIMENE': 'CC(=C)C=CCC(=C)C=C', 'B-OCIMENE': 'CC(=C)C=CCC(=C)C=C',
    'CAMPHENE': 'CC1(C)C2CCC(C2)C1=C',
    'CARENE': 'CC1=CCC2C(C1)C2(C)C', '3CARENE': 'CC1=CCC2C(C1)C2(C)C',
    '2CARENE': 'CC1=CCC2CC1C2(C)C',
    'SABINENE': 'CC(=C)C1CCC2(C)CC1C2',
    'TERPINENE': 'CC(C)C1=CC=C(C)CC1', 'ALPHATERPINENE': 'CC(C)C1=CC=C(C)CC1',
    'GAMMATERPINENE': 'CC1=CCC(=C(C)C)CC1', 'G-TERPINENE': 'CC1=CCC(=C(C)C)CC1',
    'TERPINOLENE': 'CC1=CCC(=C(C)C)CC1',
    'PCYMENE': 'Cc1ccc(C(C)C)cc1', 'PARACYMENE': 'Cc1ccc(C(C)C)cc1',
    'CYMENE': 'Cc1ccc(C(C)C)cc1',
    'PHELLANDRENE': 'CC(C)C1CCC(=C)C=C1', 'ALPHAPHELLANDRENE': 'CC(C)C1CCC(=C)C=C1',
    'BETAPHELLANDRENE': 'CC(C)=C1CCC(=C)CC1',
    # Terpene oxidation products
    'PINONALDEHYDE': 'CC1(C)C(C=O)CC(=O)C1', 'PINONAL': 'CC1(C)C(C=O)CC(=O)C1',
    'PINONICACID': 'CC1(C)C(C(=O)O)CC(=O)C1', 'PINONIC': 'CC1(C)C(C(=O)O)CC(=O)C1',
    'NORPINONICACID': 'CC1(C)C(C(=O)O)CCC1=O',
    'PINICACID': 'OC(=O)C1CC(C(=O)O)C1(C)C', 'PINIC': 'OC(=O)C1CC(C(=O)O)C1(C)C',
    'NOPINONE': 'CC1(C)C2CCC(=O)C1C2',
    'CAMPHOR': 'CC1(C)C2CCC1(C)C(=O)C2',
    'VERBENONE': 'CC1=CC(=O)C2CC1C2(C)C',
    'MYRTENOL': 'CC1=CCC2CC1C2(C)CO',
    'MYRTENAL': 'CC1=CCC2CC1C2(C)C=O',
    'TERPINEOL': 'CC(C)=C1CCC(C)(O)CC1', 'ALPHATERPINEOL': 'CC(C)=C1CCC(C)(O)CC1',
    'LINALOOL': 'CC(=C)CCCC(C)(O)C=C',
    'GERANIOL': 'CC(=C)CCCC(C)=CCO',
    'CITRONELLOL': 'CC(=C)CCCC(C)CCO',
    'CITRAL': 'CC(=C)CCCC(C)=CC=O', 'GERANIAL': 'CC(=C)CCCC(C)=CC=O',
    'CITRONELLAL': 'CC(=C)CCCC(C)CC=O',
    'MENTHOL': 'CC(C)C1CCC(C)CC1O',
    'MENTHONE': 'CC(C)C1CCC(C)CC1=O',
    'CARVONE': 'CC(=C)C1CC=C(C)C(=O)C1',
    'CARVEOL': 'CC(=C)C1CC=C(C)C(O)C1',
    'EUCALYPTOL': 'CC1(C)OC2CCC1(C)CC2', '1,8CINEOLE': 'CC1(C)OC2CCC1(C)CC2',
    'BORNEOL': 'CC1(C)C2CCC1(C)C(O)C2',
    'FENCHONE': 'CC1(C)C2CCC(C)(C2)C1=O',
    'THUJONE': 'CC(C)C12CC(=O)C(C)C1C2',
    # Other C10
    'DECANAL': 'CCCCCCCCCC=O', 'CAPRALDEHYDE': 'CCCCCCCCCC=O',
    'DECANOL': 'CCCCCCCCCCO', '1DECANOL': 'CCCCCCCCCCO',
    'DECANOICACID': 'CCCCCCCCCC(=O)O', 'CAPRICACID': 'CCCCCCCCCC(=O)O',
    'SEBACICACID': 'OC(=O)CCCCCCCCC(=O)O', 'SEBACIC': 'OC(=O)CCCCCCCCC(=O)O',
    'NAPHTHALENE': 'c1ccc2ccccc2c1', 'NAPH': 'c1ccc2ccccc2c1',

    # =========================================================================
    # C11-C15 COMPOUNDS
    # =========================================================================
    'C11H24': 'CCCCCCCCCCC', 'UNDECANE': 'CCCCCCCCCCC', 'C11000': 'CCCCCCCCCCC',
    'C12H26': 'CCCCCCCCCCCC', 'DODECANE': 'CCCCCCCCCCCC', 'C12000': 'CCCCCCCCCCCC',
    'C13H28': 'CCCCCCCCCCCCC', 'TRIDECANE': 'CCCCCCCCCCCCC', 'C13000': 'CCCCCCCCCCCCC',
    'C14H30': 'CCCCCCCCCCCCCC', 'TETRADECANE': 'CCCCCCCCCCCCCC', 'C14000': 'CCCCCCCCCCCCCC',
    'C15H32': 'CCCCCCCCCCCCCCC', 'PENTADECANE': 'CCCCCCCCCCCCCCC', 'C15000': 'CCCCCCCCCCCCCCC',
    # Sesquiterpenes (C15H24)
    'CARYOPHYLLENE': 'CC1=CCCC(=C)C2CC(C)(C)C2CC1', 'BCARYOPHYLLENE': 'CC1=CCCC(=C)C2CC(C)(C)C2CC1',
    'BETACARYOPHYLLENE': 'CC1=CCCC(=C)C2CC(C)(C)C2CC1',
    'HUMULENE': 'CC1=CCC(C)(C)C=CCC(C)=CCC1', 'ALPHAHUMULENE': 'CC1=CCC(C)(C)C=CCC(C)=CCC1',
    'FARNESENE': 'CC(=C)CCCC(=C)CCC=C(C)C', 'ALPHAFARNESENE': 'CC(=C)CCCC(=C)CCC=C(C)C',
    'NEROLIDOL': 'CC(=C)CCCC(C)(O)CCC=C(C)C',
    'BISABOLENE': 'CC1=CCC(CC1)C(C)=CCC=C(C)C',
    'LONGIFOLENE': 'CC1(C)CC2C(CC2C(C)C)C1=C',

    # =========================================================================
    # SPECIAL SPECIES / RADICALS
    # =========================================================================
    'RO2': '[O]O', 'RO': '[O]', 'R': '[C]',
    'CH3CO3H': 'CC(=O)OO', 'PERACETIC': 'CC(=O)OO',

    # =========================================================================
    # GECKO SPECIAL CODES (additional)
    # =========================================================================
    # Ethane oxidation codes
    'C02001': 'CCO[O]',  # Ethyl peroxy
    'C02002': 'CC[O]',   # Ethoxy
    'C02003': 'CCOO',    # Ethyl hydroperoxide
    'C02004': 'CCO[N+](=O)[O-]',  # Ethyl nitrate
    # Propane codes
    'C03001': 'CCCO[O]',
    'C03002': 'CCC[O]',
    'C03003': 'CCCOO',
    # Alpha-pinene oxidation codes
    'AP01': 'CC1=CCC2CC1C2(C)C',  # Alpha-pinene
    'AP02': 'CC1(C)C2CCC(O[O])C1C2',  # Pinyl peroxy
    'AP03': 'CC1(C)C2CCC([O])C1C2',   # Pinyl alkoxy
    # Beta-pinene oxidation codes
    'BP01': 'CC1(C)C2CCC(=C)C1C2',  # Beta-pinene
    'BP02': 'CC1(C)C2CCC(CO[O])C1C2',  # Nopinyl peroxy
    # Limonene codes
    'LIM01': 'CC1=CCC(CC1)C(=C)C',
    'LIM02': 'CC1=CCC(CC1)C(C)(C)O[O]',
}


# ==============================================================================
# RDKit Validation (optional but recommended)
# ==============================================================================

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.debug("RDKit not available - SMILES validation disabled")


def validate_smiles_with_rdkit(smiles: str) -> Optional[str]:
    """
    Validate and canonicalize SMILES using RDKit.
    Returns canonical SMILES if valid, None if invalid.
    """
    if not HAS_RDKIT or not smiles:
        return smiles

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass

    return None


# ==============================================================================
# Unified SMILES Conversion Function
# ==============================================================================

def gecko_to_smiles(gecko_formula: str, species_code: str = "") -> str:
    """
    Converts GECKO-A formula notation to SMILES notation.

    This is the AUTHORITATIVE function for SMILES conversion.
    All other modules should use this function.

    Strategy (in priority order):
    1. Check species code in KNOWN_SPECIES database
    2. Check formula in KNOWN_SPECIES database
    3. Parse GECKO notation to generate SMILES
    4. Validate with RDKit if available
    5. Return fallback structure if all else fails

    Args:
        gecko_formula: The GECKO formula string (e.g., 'c1HcHcHcHcHc1CH3')
        species_code: Optional GECKO species code (e.g., 'R07000')

    Returns:
        Valid SMILES string for visualization
    """
    if not gecko_formula and not species_code:
        return ""

    # 1. Check formula in known database (highest priority if specific)
    # This handles exact matches like "CH4", "H2O", "ISOPRE"
    if gecko_formula:
        formula_clean = gecko_formula.upper().replace('-', '').replace(' ', '')
        if formula_clean in KNOWN_SPECIES:
            return KNOWN_SPECIES[formula_clean]

    # 2. Try to PARSE the formula structurally (New "Better Logic")
    # GECKO often provides structural formulas like "CH3CH(CH3)CH2CH3".
    # We should trust these over generic species codes (like C05000) which
    # might map to n-alkanes by default.
    if gecko_formula:
        # First check for aromatic notation which requires special handling
        aromatic_smiles = _parse_aromatic_gecko(gecko_formula)
        if aromatic_smiles:
            validated = validate_smiles_with_rdkit(aromatic_smiles)
            if validated:
                return validated
            return aromatic_smiles

        # Then try standard structural parsing
        parsed_smiles = _convert_gecko_to_smiles(gecko_formula)
        
        # If parsing produced a valid non-empty SMILES that looks like a structure
        if parsed_smiles and len(parsed_smiles) > 0:
            # Validate with RDKit to be sure
            validated = validate_smiles_with_rdkit(parsed_smiles)
            if validated:
                 return validated
            
            # If validation fails, do NOT return the parsed result.
            # We must fall back to the species code lookup (Step 3).
            # This prevents "boxes" caused by invalid parsed SMILES.

    # 3. Check species code in known database (Fallback)
    # Only if formula parsing failed (e.g. formula was "C5H12" or empty),
    # do we rely on the static code mapping.
    if species_code:
        code_upper = species_code.upper().strip()

        # Handle 'A' prefix (Aerosol species) which often match 'H' species or base species
        # e.g. AH05000 -> H05000
        possible_variants = [
            code_upper,
            code_upper.lstrip('G'),   # Remove G prefix
            code_upper.replace('-', '').replace('_', ''),
            code_upper.replace('ALPHA', 'A').replace('BETA', 'B')
        ]
        
        # Add fallback for Aerosol species (convert AH05000 -> H05000)
        if code_upper.startswith('A') and len(code_upper) > 1:
             possible_variants.append(code_upper[1:]) # Strip A
             if code_upper[1] == 'H':
                 possible_variants.append(code_upper[2:]) # Strip AH

        for variant in possible_variants:
            if variant in KNOWN_SPECIES:
                return KNOWN_SPECIES[variant]
                
    # 3b. (Aromatic handling moved to Step 2)
    
    # 3c. Final check: If it's an 'A' code and we still have no formula,
    # allows the system to query the database/dictionary in the caller for the stripped code.
    # (This logic is usually handled by the caller, but if we are here, we verify locally).
    
    # 4. Handle simple inorganic molecules (Final check)
    if gecko_formula in KNOWN_SPECIES:
        return KNOWN_SPECIES[gecko_formula]

    # 5. (Removed - functionality moved to Step 2)
    
    # 6. Final fallback
    # If we reached here, parsing failed and code lookup failed.
    # We already tried parsing in Step 2, but _convert_gecko_to_smiles might return ""
    # for things it doesn't understand.
    
    fallback = _generate_fallback_smiles(gecko_formula)
    logger.info(f"Using fallback SMILES for '{gecko_formula}' (code: {species_code}): {fallback}")
    return fallback


def _parse_aromatic_gecko(formula: str) -> Optional[str]:
    """
    Parse GECKO aromatic ring notation.

    GECKO uses patterns like:
    - c1HcHcHcHcHc1H (benzene)
    - c1(CH3)cHcHcHcHc1H (toluene with substituent at position 1)
    - c1HcHcHcHcHc1CH3 (toluene with substituent after ring)
    """
    if not formula or 'c' not in formula.lower():
        return None

    # Pattern 1: c1(substituent)cHcHcHcHc1H
    pattern1 = re.search(r'c(\d)\(([^)]+)\)(?:cH){4}c\1H?', formula)
    if pattern1:
        ring_num = pattern1.group(1)
        substituent_gecko = pattern1.group(2)
        substituent = _convert_substituent(substituent_gecko)
        if substituent:
            return f"{substituent}c1ccccc1"
        return "c1ccccc1"

    # Pattern 2: c1HcHcHcHcHc1(substituent) or c1HcHcHcHcHc1substituent
    pattern2 = re.search(r'c(\d)H(?:cH){4}c\1(.*)$', formula)
    if pattern2:
        rest = pattern2.group(2)
        # Check for parenthesized substituent
        if rest.startswith('(') and ')' in rest:
            sub_end = rest.index(')')
            substituent_gecko = rest[1:sub_end]
        else:
            substituent_gecko = rest
        substituent = _convert_substituent(substituent_gecko)
        if substituent:
            return f"{substituent}c1ccccc1"
        return "c1ccccc1"

    # Pattern 3: Simple benzene c1HcHcHcHcHc1 or c1cHcHcHcHc1H
    if re.search(r'c\d[Hc]{8,}c\d', formula):
        return "c1ccccc1"

    return None


def _convert_substituent(rest: str) -> str:
    """Convert GECKO substituent notation to SMILES.

    Handles substituent patterns for aromatic compounds like ethylbenzene derivatives.

    Examples:
        CH2CH3 -> CC (ethyl)
        CH(OO.)CH3 -> CC(O[O]) (1-phenylethyl peroxy radical)
        CH(O.)CH3 -> CC([O]) (1-phenylethyl alkoxy radical)
        CH(OH)CH3 -> CC(O) (1-phenylethanol)
        CH(ONO2)CH3 -> CC(O[N+](=O)[O-]) (1-phenylethyl nitrate)
        CH(OOH)CH3 -> CC(OO) (1-phenylethyl hydroperoxide)
        COCH3 -> C(=O)C (acetyl - acetophenone)
        CH3 -> C (methyl)
    """
    if not rest:
        return ""

    rest = rest.strip()

    # =========================================================================
    # PATTERN 1: CH(functional_group)CH3 - substituted ethyl groups
    # This is CRITICAL for ethylbenzene derivatives (R08000 series)
    # The pattern is: CH-carbon with functional group attached, then terminal CH3
    # =========================================================================

    # CH(OO.)CH3 -> peroxy radical on alpha carbon of ethyl group
    if rest.startswith('CH(OO.)CH3') or rest.startswith('CH(OO)CH3'):
        return 'CC(O[O])'  # ethyl with peroxy on alpha carbon

    # CH(O.)CH3 -> alkoxy radical on alpha carbon of ethyl group
    if rest.startswith('CH(O.)CH3'):
        return 'CC([O])'  # ethyl with alkoxy radical on alpha carbon

    # CH(OH)CH3 -> alcohol on alpha carbon of ethyl group (1-phenylethanol)
    if rest.startswith('CH(OH)CH3'):
        return 'CC(O)'  # ethyl with hydroxyl on alpha carbon

    # CH(ONO2)CH3 -> nitrate on alpha carbon of ethyl group
    if rest.startswith('CH(ONO2)CH3'):
        return 'CC(O[N+](=O)[O-])'  # ethyl with nitrate on alpha carbon

    # CH(OOH)CH3 -> hydroperoxide on alpha carbon of ethyl group
    if rest.startswith('CH(OOH)CH3'):
        return 'CC(OO)'  # ethyl with hydroperoxide on alpha carbon

    # =========================================================================
    # PATTERN 2: CO patterns - ketones/carbonyls
    # =========================================================================

    # COCH3 -> acetyl group (acetophenone when on benzene)
    if rest.startswith('COCH3') or rest.startswith('(CO)CH3'):
        return 'C(=O)C'  # acetyl group

    # CO alone or (CO) -> aldehyde/carbonyl
    if rest == 'CO' or rest == '(CO)':
        return 'C=O'

    # =========================================================================
    # PATTERN 3: Simple alkyl chains (check longer ones first!)
    # =========================================================================

    # Ethyl group: CH2CH3, C2H5
    if rest.startswith('CH2CH3'):
        return 'CC'
    if rest.startswith('C2H5') or rest.startswith('C2H4'):
        return 'CC'

    # Propyl group: CH2CH2CH3, C3H7
    if rest.startswith('CH2CH2CH3'):
        return 'CCC'
    if rest.startswith('C3H7') or rest.startswith('C3H6'):
        return 'CCC'

    # Isopropyl: CH(CH3)2 or (CH3)2CH
    if rest.startswith('CH(CH3)2') or rest.startswith('(CH3)2CH'):
        return 'C(C)C'

    # Butyl group: C4H9
    if rest.startswith('C4H9') or rest.startswith('C4H8'):
        return 'CCCC'

    # =========================================================================
    # PATTERN 4: Simple terminal groups (shortest patterns last)
    # =========================================================================

    # Methyl group: CH3
    if rest.startswith('CH3'):
        return 'C'

    # Aldehyde: CHO
    if rest.startswith('CHO'):
        return 'C=O'

    # Methylene: CH2 (when not followed by more carbon chain)
    if rest.startswith('CH2') and not (len(rest) > 3 and rest[3:].startswith('CH')):
        return 'C'

    # =========================================================================
    # PATTERN 5: Handle CH(X) patterns where X is a functional group
    # This catches remaining patterns not covered above
    # =========================================================================

    # Generic pattern: CH(functional)remainder
    ch_func_match = re.match(r'^CH\(([^)]+)\)(.*)$', rest)
    if ch_func_match:
        func_group = ch_func_match.group(1)
        remainder = ch_func_match.group(2)

        # Convert the functional group
        func_smiles = ''
        if func_group in ('OO.', 'OO'):
            func_smiles = 'O[O]'
        elif func_group == 'O.':
            func_smiles = '[O]'
        elif func_group == 'OH':
            func_smiles = 'O'
        elif func_group == 'ONO2':
            func_smiles = 'O[N+](=O)[O-]'
        elif func_group == 'OOH':
            func_smiles = 'OO'
        elif func_group == 'CO':
            func_smiles = '=O'

        # Convert the remainder (usually CH3 for ethyl derivatives)
        remainder_smiles = ''
        if remainder:
            if remainder.startswith('CH3'):
                remainder_smiles = 'C'
            elif remainder.startswith('CH2CH3'):
                remainder_smiles = 'CC'
            elif remainder.startswith('CH2'):
                remainder_smiles = 'C'

        if func_smiles:
            return f'{remainder_smiles}C({func_smiles})'

    # Single CH (not CHO)
    if rest.startswith('CH') and not rest.startswith('CHO'):
        return 'C'

    # If nothing matched, return empty (fallback to benzene)
    return ""


def _convert_gecko_to_smiles(gecko_formula: str) -> str:
    """
    Convert GECKO structural notation to SMILES.
    Handles non-aromatic compounds.
    """
    s = gecko_formula

    # Remove radical markers (dots after atoms in parentheses)
    s = re.sub(r'(\w+)\.\)', r'\1)', s)
    s = s.replace('.)', ')')

    # Convert GECKO structural notations step by step
    # IMPORTANT: Preserve outer parentheses for branching structure

    # Handle nitrate groups: (ONO2) -> (O[N+](=O)[O-])
    s = s.replace('(ONO2)', '(O[N+](=O)[O-])')
    s = s.replace('ONO2', 'O[N+](=O)[O-]') # If not in parens

    # Handle hydroperoxides: (OOH) -> (OO)
    s = s.replace('(OOH)', '(OO)')
    s = s.replace('OOH', 'OO')

    # Handle peroxy radical: (OO.) -> (O[O])
    s = s.replace('(OO.)', '(O[O])')
    s = s.replace('(OO)', '(OO)')

    # Handle alkoxy radical: (O.) -> ([O])
    s = s.replace('(O.)', '([O])')

    # Handle hydroxyl: (OH) -> (O)
    s = s.replace('(OH)', '(O)')

    # Handle carbonyl: (CO) -> (=O)
    s = s.replace('(CO)', '(=O)')

    # Aldehydes: CHO at end -> C=O
    s = re.sub(r'CHO$', 'C=O', s)
    s = re.sub(r'CHO([^A-Za-z])', r'C=O\1', s)

    # Convert explicit hydrogens: CH3, CH2, CH -> C
    s = s.replace('CH3', 'C')
    s = s.replace('CH2', 'C')
    # Be careful with CH - don't replace CHO which is already handled
    s = re.sub(r'CH(?![OA-Z])', 'C', s)

    # Handle double bonds: Cd=Cd or CdH=CdH2 patterns
    s = s.replace('Cd=Cd', 'C=C')
    s = s.replace('CdH2', 'C')
    s = s.replace('CdH', 'C')
    s = s.replace('Cd', 'C')

    # Clean up invalid patterns
    s = s.replace('===', '=')
    s = s.replace('==', '=')

    # Remove any remaining dots
    s = re.sub(r'\.(?=[^a-z]|$)', '', s)

    # Check if result looks invalid (unclosed brackets, etc.)
    open_parens = s.count('(')
    close_parens = s.count(')')
    if open_parens != close_parens:
        # Unbalanced parentheses - fall back to simple representation
        return _generate_fallback_smiles(gecko_formula)

    return s


def _generate_fallback_smiles(formula: str) -> str:
    """
    Generate a representative fallback SMILES when parsing fails.
    """
    # Count approximate carbons
    n_carbon = len(re.findall(r'C[dH123]*', formula))
    if n_carbon == 0:
        match = re.search(r'C(\d+)', formula)
        if match:
            n_carbon = int(match.group(1))

    n_carbon = max(1, min(n_carbon, 15))

    # Detect features
    has_double = 'Cd' in formula or '=C' in formula
    has_ring = bool(re.search(r'[CO]\d.*[CO]\d', formula))
    has_aromatic = 'cH' in formula.lower()

    has_oh = '(OH)' in formula
    has_ooh = '(OOH)' in formula
    has_ono2 = '(ONO2)' in formula
    has_peroxy = '(OO.)' in formula or 'OO.' in formula
    has_alkoxy = '(O.)' in formula and not has_peroxy
    has_carbonyl = '(CO)' in formula or '-CO-' in formula
    has_aldehyde = 'CHO' in formula

    # Build base structure
    if has_aromatic and n_carbon >= 6:
        base = 'c1ccccc1'
        if n_carbon > 6:
            base = 'C' * (n_carbon - 6) + base
    elif has_ring and n_carbon >= 4:
        ring_size = min(6, n_carbon)
        base = 'C1' + 'C' * (ring_size - 2) + 'C1'
        if n_carbon > ring_size:
            base += 'C' * (n_carbon - ring_size)
    elif has_double and n_carbon >= 2:
        if n_carbon == 5:
            base = 'CC(=C)C=C'  # Isoprene-like
        elif n_carbon >= 4:
            half = n_carbon // 2
            base = 'C' * half + '=C' + 'C' * (n_carbon - half - 1)
        else:
            base = 'C=C'
    else:
        base = 'C' * n_carbon

    # Add functional groups
    suffix = ''
    if has_oh:
        suffix += 'O'
    if has_ooh:
        suffix += 'OO'
    if has_ono2:
        suffix += 'O[N+](=O)[O-]'
    if has_peroxy:
        suffix += 'O[O]'
    if has_alkoxy:
        suffix += '[O]'

    if has_carbonyl or has_aldehyde:
        if len(base) >= 3:
            mid = len(base) // 2
            base = base[:mid] + '(=O)' + base[mid:]
        else:
            suffix += 'C=O'

    return base + suffix


# ==============================================================================
# Reaction Tree Parsing
# ==============================================================================

def parse_reaction_tree(output_dir: str, voc_name: str, max_depth: int = 3, max_nodes: int = 50) -> Dict:
    """
    Parses the mechanism file to build a reaction tree starting from the VOC.
    Returns a JSON-serializable dictionary representing the graph.

    Args:
        output_dir: Directory containing GECKO-A output files
        voc_name: Name of the VOC (e.g., 'toluene', 'isoprene')
        max_depth: Maximum depth of the tree traversal
        max_nodes: Maximum number of nodes to include

    Returns:
        Dictionary with 'nodes' and 'edges' lists for Cytoscape visualization
    """
    # 1. Find the mechanism file
    mech_file = os.path.join(output_dir, "reactions.txt")
    if not os.path.exists(mech_file):
        import glob
        mecs = glob.glob(os.path.join(output_dir, "*.mec"))
        if mecs:
            mech_file = mecs[0]
        else:
            logger.warning("No mechanism file found for tree generation.")
            return {"nodes": [], "edges": []}

    # 2. Read and parse the dictionary file for species info
    code_map = _read_dictionary(output_dir)

    # 3. Find the root node
    root_node = _find_root_node(output_dir, voc_name, code_map)
    logger.info(f"Using root node: {root_node} for VOC: {voc_name}")

    # 4. Parse reactions and build adjacency list
    adj_list = _parse_reactions(mech_file, code_map)

    # 5. Handle case where root is not found
    if root_node not in adj_list:
        root_node = _find_alternative_root(adj_list, root_node)
        if root_node not in adj_list:
            logger.warning(f"Root node {root_node} has no edges. Tree will be empty.")
            return {"nodes": [], "edges": []}

    # 6. Build tree using BFS
    nodes, edges = _build_tree_bfs(adj_list, root_node, max_depth, max_nodes)

    # 7. Format for Cytoscape with SMILES
    node_list = _format_nodes(nodes, code_map, voc_name)
    edge_list = _format_edges(edges, nodes, code_map)

    logger.info(f"Reaction tree generated: {len(node_list)} nodes, {len(edge_list)} edges")
    return {"nodes": node_list, "edges": edge_list}


def _read_dictionary(output_dir: str) -> Dict[str, str]:
    """Read species codes and formulas from dictionary.out"""
    code_map = {}
    dict_path = os.path.join(output_dir, "dictionary.out")

    if os.path.exists(dict_path):
        with open(dict_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    code_map[parts[0]] = parts[1]

    return code_map


def _find_root_node(output_dir: str, voc_name: str, code_map: Dict[str, str]) -> str:
    """
    Find the root VOC node using multiple strategies.
    Priority order:
    1. Read from listprimary.dat (authoritative)
    2. Exact match in code_map
    3. Fuzzy match on code
    4. Heuristic based on carbon count
    """
    # 1. Read from listprimary.dat (most authoritative)
    listprimary_path = os.path.join(output_dir, "listprimary.dat")
    if os.path.exists(listprimary_path):
        try:
            with open(listprimary_path, 'r') as f:
                for line in f:
                    code = line.strip().split()[0] if line.strip() else None
                    if code and code in code_map:
                        logger.info(f"Root VOC from listprimary.dat: {code}")
                        return code
        except Exception as e:
            logger.debug(f"Could not read listprimary.dat: {e}")

    # 2. Exact match on formula/name
    for code, formula in code_map.items():
        if formula.lower() == voc_name.lower():
            return code

    # 3. Fuzzy match on code
    v_upper = voc_name.upper().replace('-', '').replace('_', '')
    for code in code_map:
        c_upper = code.upper()
        if (v_upper.startswith(c_upper) and len(c_upper) >= 4) or \
           (c_upper.startswith(v_upper) and len(v_upper) >= 4):
            return code

    # 4. Special case mappings
    special_mappings = {
        'ISOPRENE': 'ISOPRE',
        'TOLUENE': 'R07000',
        'ETHANE': 'C02000',
        'PROPANE': 'C03000',
        'BUTANE': 'C04000',
        'PENTANE': 'C05000',
        'NPENTANE': 'C05000',
        'ISOPENTANE': 'IPENTA',  # 2-methylbutane
        '2METHYLBUTANE': 'IPENTA',
        'NEOPENTANE': 'NEOPEN',  # 2,2-dimethylpropane
        '22DIMETHYLPROPANE': 'NEOPEN',
        'ISOBUTANE': 'IBUTAN',  # 2-methylpropane
        '2METHYLPROPANE': 'IBUTAN',
        'HEXANE': 'C06000',
        'NHEXANE': 'C06000',
        'CYCLOHEXANE': 'CYCHEX',
        'CYCLOPENTANE': 'CYCPEN',
        'ALPHAPINENE': 'APINEN',
        'ALPHA-PINENE': 'APINEN',
        'BETAPINENE': 'BPINEN',
        'BETA-PINENE': 'BPINEN',
        'LIMONENE': 'LIMONE',
        'MYRCENE': 'MYRCEN',
        'OCIMENE': 'OCIMEN',
    }

    voc_key = voc_name.upper().replace('-', '').replace('_', '')
    if voc_key in special_mappings:
        mapped = special_mappings[voc_key]
        if mapped in code_map:
            return mapped

    # 5. Heuristic: First organic species with C >= 2
    dict_path = os.path.join(output_dir, "dictionary.out")
    if os.path.exists(dict_path):
        try:
            with open(dict_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) > 13:
                        code = parts[0]
                        try:
                            nC = int(parts[3])
                            if nC >= 2:
                                logger.info(f"Heuristic root node found: {code} (nC={nC})")
                                return code
                        except (ValueError, IndexError):
                            pass
        except Exception as e:
            logger.debug(f"Heuristic search failed: {e}")

    return voc_name


def _find_alternative_root(adj_list: Dict, original_root: str) -> str:
    """Find an alternative root when the original is not in the adjacency list."""
    if not adj_list:
        return original_root

    # Find source nodes (nodes that appear only as reactants, not products)
    all_products = set()
    for src, dests in adj_list.items():
        for d in dests:
            all_products.add(d['to'])

    potential_roots = [n for n in adj_list.keys() if n not in all_products]

    # Filter out inorganic species
    inorganics = {"OH", "HO2", "NO", "NO2", "O3", "HV", "M", "O2", "H2O", "CO", "CO2"}
    potential_roots = [n for n in potential_roots if n not in inorganics]

    if potential_roots:
        # Pick the one with the shortest name
        root = min(potential_roots, key=len)
        logger.info(f"Found alternative root: {root}")
        return root

    # Fallback: node with most edges
    if adj_list:
        root = max(adj_list.keys(), key=lambda k: len(adj_list[k]))
        logger.info(f"Fallback to most connected node: {root}")
        return root

    return original_root


def _parse_reactions(mech_file: str, code_map: Dict[str, str]) -> Dict:
    """Parse reaction file and build adjacency list."""
    adj_list = {}

    with open(mech_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        # Skip comments
        if line.strip().startswith("!"):
            continue

        clean_line = line.split("!")[0].strip()

        # Determine separator
        sep = None
        if "=>" in clean_line:
            sep = "=>"
        elif "=" in clean_line:
            sep = "="

        if not sep:
            continue

        parts = clean_line.split(sep)
        if len(parts) != 2:
            continue

        lhs = parts[0].strip()
        rhs = parts[1].strip()

        # Parse reactants and products
        reactants = _parse_species_list(lhs, code_map)
        products = _parse_species_list(rhs, code_map, parse_yields=True)

        # Add edges
        for reactant in reactants:
            if reactant not in adj_list:
                adj_list[reactant] = []
            for product, yield_val in products:
                adj_list[reactant].append({'to': product, 'yield': yield_val})

    return adj_list


def _parse_species_list(species_str: str, code_map: Dict, parse_yields: bool = False) -> List:
    """Parse a species string (LHS or RHS of reaction)."""
    species = []

    for token in species_str.split("+"):
        token = token.strip()
        if not token:
            continue

        coeff = 1.0
        name = token

        # Check for leading coefficient
        match = re.match(r'^([\d\.]+)\s+(.*)$', token)
        if match:
            try:
                coeff = float(match.group(1))
                name = match.group(2)
            except ValueError:
                pass

        # Extract species name
        parts = name.split()
        if parts:
            name = parts[0]

        # Normalize G-prefix
        if name.startswith("G") and len(name) > 1:
            candidate = name[1:]
            if candidate in code_map or candidate in KNOWN_SPECIES:
                name = candidate

        # Filter non-species
        non_species = {"HV", "FALLOFF", "EXTRA", "NOTHING", "TBODY", "OXYGEN", "M"}
        if name in non_species or not re.search(r'[a-zA-Z]', name):
            continue

        if parse_yields:
            species.append((name, coeff))
        else:
            species.append(name)

    return species


def _build_tree_bfs(adj_list: Dict, root: str, max_depth: int, max_nodes: int) -> Tuple[set, List]:
    """Build reaction tree using BFS traversal."""
    nodes = {root}
    edges = []
    visited = {root}
    queue = [(root, 0)]

    # Track seen edges to avoid duplicates (use dict to keep highest yield)
    seen_edges: Dict[Tuple[str, str], float] = {}

    # Inorganic species to filter out
    inorganics = {"OH", "HO2", "NO", "NO2", "O3", "CO", "CO2", "H2O", "O2", "SO2", "SULF", "H2"}

    while queue and len(nodes) < max_nodes:
        current, depth = queue.pop(0)

        if depth >= max_depth:
            continue

        if current not in adj_list:
            continue

        for child_obj in adj_list[current]:
            child = child_obj['to']
            yield_val = child_obj['yield']

            # Filter inorganics
            if child in inorganics:
                continue

            # Filter low yields
            if yield_val < 0.05:
                continue

            # Deduplicate edges - keep the one with highest yield
            edge_key = (current, child)
            if edge_key not in seen_edges or yield_val > seen_edges[edge_key]:
                seen_edges[edge_key] = yield_val

            nodes.add(child)

            if child not in visited:
                visited.add(child)
                queue.append((child, depth + 1))

    # Convert seen_edges dict to list
    for (from_node, to_node), yield_val in seen_edges.items():
        edges.append({"from": from_node, "to": to_node, "yield": yield_val})

    return nodes, edges


def identify_functional_groups(smiles: str, formula: str = "") -> List[str]:
    """Identify functional groups from SMILES and formula."""
    groups = []
    
    if not smiles:
        return groups

    # Check for Nitrates
    if 'N+' in smiles and 'O-' in smiles: # Standard RDKit representation
        groups.append("Nitrate")
    elif 'ONO2' in formula:
        groups.append("Nitrate")

    # Check for PANs
    if 'C(=O)OO[N+]' in smiles or 'C(=O)OON' in smiles:
        groups.append("PAN-like")
    elif 'C(=O)OO' in smiles:
         # Peroxy acid or similar
         if 'Nitrate' not in groups: # Don't double count PAN
             groups.append("Peroxy Acid")

    # Check for Hydroperoxides
    if 'OO' in smiles and 'OO[N+]' not in smiles and 'C(=O)OO' not in smiles:
        groups.append("Hydroperoxide")
    elif 'OOH' in formula:
        groups.append("Hydroperoxide")
        
    # Check for Alcohols (OH)
    # This is tricky with SMILES as 'O' is ether, alcohol, etc.
    # We look for explicit [OH] or standard O patterns if simple
    if '[OH]' in smiles:
         groups.append("Alcohol")
    elif 'OH' in formula and 'COOH' not in formula and 'OOH' not in formula:
         groups.append("Alcohol")
         
    # Check for Carbonyls
    if '=O' in smiles:
        if 'C=O' in smiles or 'C(=O)' in smiles:
             # Distinguish Aldehyde vs Ketone vs Acid is hard without graph analysis
             # But we can check formula for hints
             if 'CHO' in formula:
                 groups.append("Aldehyde")
             elif 'COOH' in formula:
                 groups.append("Carboxylic Acid")
             elif 'CO' in formula and 'COOH' not in formula:
                 groups.append("Ketone")
             else:
                 groups.append("Carbonyl")
    
    # Check for Epoxides
    if 'C1OC1' in smiles or 'C1OC1' in formula or 'EPOX' in formula:
        groups.append("Epoxide")

    # Unique filtering
    return list(set(groups))


def _format_nodes(nodes: set, code_map: Dict, voc_name: str = "") -> List[Dict]:
    """Format nodes for Cytoscape visualization."""
    node_list = []

    is_isopentane_job = False
    if voc_name:
        voc_clean = voc_name.lower().strip().replace('-', '').replace('_', '')
        if voc_clean in ['isopentane', '2methylbutane', 'ipentane']:
            is_isopentane_job = True

    for n in nodes:
        raw_formula = code_map.get(n, n)
        smiles = gecko_to_smiles(raw_formula, species_code=n)

        # --- FIX FOR ISOPENTANE/C05000 AMBIGUITY ---
        if is_isopentane_job:
            if (n in {'C05000', 'G05000', '205000', 'GH05000', 'H05000'}) and (smiles in {'CCCCC', 'CCC(C)C'} or not smiles):
                smiles = 'CC(C)CC'
        # -------------------------------------------

        label = raw_formula[:17] + "..." if len(raw_formula) > 20 else raw_formula
        
        # Identify functional groups
        fgs = identify_functional_groups(smiles, raw_formula)

        node_list.append({
            "id": n,
            "label": label,
            "smiles": smiles,
            "raw_formula": raw_formula,
            "code": n,
            "display_code": n,
            "title": f"{n}: {raw_formula}",
            "functional_groups": fgs  # Added field
        })

    return node_list


def _format_edges(edges: List[Dict], nodes: set, code_map: Dict) -> List[Dict]:
    """Format edges for Cytoscape visualization."""
    enhanced_edges = []

    for edge in edges:
        # Only include edges where both endpoints are in nodes
        if edge['from'] not in nodes or edge['to'] not in nodes:
            continue

        from_formula = code_map.get(edge['from'], edge['from'])
        to_formula = code_map.get(edge['to'], edge['to'])

        enhanced_edges.append({
            "from": edge['from'],
            "to": edge['to'],
            "yield": edge['yield'],
            "reaction": f"{from_formula} → {to_formula}",
            "label": f"{edge['yield']:.2f}" if edge['yield'] != 1.0 else ""
        })

    return enhanced_edges


if __name__ == "__main__":
    # Test the SMILES conversion
    test_cases = [
        ("c1HcHcHcHcHc1H", ""),  # Benzene
        ("c1(CH3)cHcHcHcHc1H", ""),  # Toluene
        ("", "R07000"),  # Toluene by code
        ("CH3Cd(=CdH2)CdH=CdH2", "ISOPRE"),  # Isoprene
        ("CC1=CCC2CC1C2(C)C", "APINEN"),  # Alpha-pinene
    ]

    for formula, code in test_cases:
        result = gecko_to_smiles(formula, code)
        print(f"Formula: {formula or 'N/A'}, Code: {code or 'N/A'} -> SMILES: {result}")
