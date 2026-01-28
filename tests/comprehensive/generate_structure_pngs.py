#!/usr/bin/env python3
"""
Generate PNG images for 50 compounds and create a validation report.

This script:
1. Generates 2D structure images from SMILES
2. Generates images from GECKO formula -> SMILES conversion
3. Creates a comparison report
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from PIL import Image
import json

from gecko_web.chemdata.compound_database import COMPOUNDS
from gecko_web.reaction_tree import gecko_to_smiles, KNOWN_SPECIES

# Output directory
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'structure_images')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 50 compounds to validate
VALIDATION_COMPOUNDS = [
    # Terpene oxidation products (CRITICAL - cyclobutane)
    'pinic_acid',
    'pinonic_acid',
    'pinonaldehyde',
    'nopinone',

    # Monoterpenes (bicyclic)
    'alpha_pinene',
    'beta_pinene',
    'camphene',
    'sabinene',
    '3_carene',

    # Monoterpenes (monocyclic)
    'limonene',
    'alpha_terpinene',
    'gamma_terpinene',
    'terpinolene',

    # Monoterpenes (acyclic)
    'myrcene',
    'ocimene',

    # Oxygenated terpenes
    'linalool',
    'alpha_terpineol',
    'nerolidol',
    '1_8_cineole',

    # Sesquiterpenes
    'beta_caryophyllene',
    'alpha_humulene',

    # Simple alkanes
    'methane',
    'ethane',
    'propane',
    'isobutane',
    'isopentane',
    'neopentane',
    'cyclohexane',

    # Alkenes
    'ethene',
    'propene',
    'isoprene',
    'isobutene',

    # Aromatics
    'benzene',
    'toluene',
    'ethylbenzene',
    'styrene',
    'naphthalene',

    # Oxygenated VOCs
    'methanol',
    'ethanol',
    'acetone',
    'formaldehyde',
    'acetaldehyde',
    'acetic_acid',
    'formic_acid',

    # Esters
    'methyl_formate',
    'ethyl_acetate',
    'butyl_acetate',

    # Ethers
    'dimethyl_ether',
    'diethyl_ether',
    'methyl_tert_butyl_ether',
]

def generate_structure_image(smiles: str, filename: str, title: str = "") -> bool:
    """Generate a PNG image from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Generate 2D coordinates
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)

        # Draw molecule
        img = Draw.MolToImage(mol, size=(400, 300), legend=title)
        img.save(filename)
        return True
    except Exception as e:
        print(f"Error generating image for {title}: {e}")
        return False

def get_ring_info(smiles: str) -> dict:
    """Get ring information from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"valid": False, "rings": [], "num_rings": 0}

    ring_info = mol.GetRingInfo()
    ring_sizes = sorted([len(r) for r in ring_info.AtomRings()])

    return {
        "valid": True,
        "rings": ring_sizes,
        "num_rings": len(ring_sizes),
        "carbons": sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C'),
        "mw": round(Descriptors.MolWt(mol), 2),
        "formula": rdMolDescriptors.CalcMolFormula(mol)
    }

def main():
    report = {
        "total_compounds": len(VALIDATION_COMPOUNDS),
        "successful": 0,
        "failed": 0,
        "compounds": []
    }

    print("=" * 70)
    print("GENERATING STRUCTURE IMAGES FOR 50 COMPOUNDS")
    print("=" * 70)

    for i, compound_name in enumerate(VALIDATION_COMPOUNDS, 1):
        print(f"\n[{i:2d}/50] Processing: {compound_name}")

        compound = COMPOUNDS.get(compound_name)
        if not compound:
            print(f"  ERROR: Compound '{compound_name}' not found in database!")
            report["failed"] += 1
            report["compounds"].append({
                "name": compound_name,
                "status": "NOT_FOUND",
                "error": "Not in compound database"
            })
            continue

        # Get compound info
        smiles = compound.smiles
        gecko_formula = compound.gecko_formula

        print(f"  SMILES: {smiles}")
        print(f"  GECKO:  {gecko_formula}")

        # Analyze structure
        info = get_ring_info(smiles)
        print(f"  Formula: {info.get('formula', 'N/A')}")
        print(f"  MW: {info.get('mw', 'N/A')}")
        print(f"  Rings: {info.get('rings', [])}")
        print(f"  Carbons: {info.get('carbons', 'N/A')}")

        # Generate image from database SMILES
        img_path = os.path.join(OUTPUT_DIR, f"{compound_name}.png")
        title = f"{compound_name}\n{info.get('formula', '')} MW={info.get('mw', '')}"
        success = generate_structure_image(smiles, img_path, compound_name)

        if success:
            print(f"  Image saved: {img_path}")
            report["successful"] += 1
        else:
            print(f"  ERROR: Failed to generate image")
            report["failed"] += 1

        # Also test gecko_to_smiles conversion
        gecko_smiles = gecko_to_smiles(gecko_formula, "")
        gecko_info = get_ring_info(gecko_smiles)

        # Check if conversion produces same ring structure
        rings_match = info.get('rings', []) == gecko_info.get('rings', [])

        report["compounds"].append({
            "name": compound_name,
            "status": "OK" if success else "FAILED",
            "smiles": smiles,
            "gecko_formula": gecko_formula,
            "gecko_to_smiles": gecko_smiles,
            "formula": info.get('formula'),
            "mw": info.get('mw'),
            "rings": info.get('rings'),
            "carbons": info.get('carbons'),
            "gecko_rings": gecko_info.get('rings'),
            "rings_match": rings_match,
            "image_path": img_path if success else None
        })

    # Special check for critical terpene oxidation products
    print("\n" + "=" * 70)
    print("CRITICAL CHECK: TERPENE OXIDATION PRODUCTS")
    print("=" * 70)

    critical = ['pinic_acid', 'pinonic_acid', 'pinonaldehyde']
    for name in critical:
        c = next((x for x in report["compounds"] if x["name"] == name), None)
        if c:
            has_4_ring = 4 in c.get('rings', [])
            has_5_ring = 5 in c.get('rings', [])
            has_6_ring = 6 in c.get('rings', [])

            status = "CORRECT (cyclobutane)" if has_4_ring and not has_5_ring and not has_6_ring else "WRONG!"
            print(f"\n{name}:")
            print(f"  Rings: {c.get('rings')}")
            print(f"  Status: {status}")
            print(f"  Formula: {c.get('formula')}")

    # Save report
    report_path = os.path.join(OUTPUT_DIR, "validation_report.json")
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved to: {report_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total compounds: {report['total_compounds']}")
    print(f"Successful: {report['successful']}")
    print(f"Failed: {report['failed']}")

    # List all 50 compounds
    print("\n" + "=" * 70)
    print("ALL 50 COMPOUNDS TESTED:")
    print("=" * 70)
    for i, c in enumerate(report["compounds"], 1):
        status = "✓" if c["status"] == "OK" else "✗"
        rings = c.get('rings', [])
        print(f"{i:2d}. {status} {c['name']:25s} rings={rings} formula={c.get('formula', 'N/A')}")

    return report

if __name__ == '__main__':
    main()
