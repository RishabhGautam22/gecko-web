#!/usr/bin/env python3
"""
Generate labeled PNG images for ALL 130 compounds in the dropdown list.
Validates each structure against PubChem/database specifications.
"""

import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, rdMolDescriptors
from PIL import Image, ImageDraw, ImageFont
import io

from gecko_web.chemdata.compound_database import COMPOUNDS

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'all_130_compounds')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def generate_labeled_structure(compound_name: str, smiles: str, formula: str,
                                output_path: str, size=(250, 200)) -> bool:
    """Generate a structure image with compound name and formula labels."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Calculate dimensions
        label_height = 40
        mol_height = size[1] - label_height

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Draw molecule
        drawer = Draw.MolDraw2DCairo(size[0], mol_height)
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        opts.bondLineWidth = 2.5
        opts.padding = 0.08
        opts.minFontSize = 12

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Load as PIL Image
        mol_img = Image.open(io.BytesIO(drawer.GetDrawingText()))

        # Create combined image with labels
        combined = Image.new('RGB', size, 'white')
        combined.paste(mol_img, (0, 0))

        # Add labels
        draw = ImageDraw.Draw(combined)

        # Try to load a font
        font = None
        font_small = None
        font_paths = [
            '/System/Library/Fonts/Helvetica.ttc',
            '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
            'C:/Windows/Fonts/arial.ttf',
        ]
        for fp in font_paths:
            if os.path.exists(fp):
                try:
                    font = ImageFont.truetype(fp, 11)
                    font_small = ImageFont.truetype(fp, 10)
                    break
                except Exception:
                    pass

        if font is None:
            font = ImageFont.load_default()
            font_small = font

        y_pos = mol_height + 3

        # Draw compound name
        display_name = compound_name.replace('_', ' ')
        if len(display_name) > 25:
            display_name = display_name[:23] + '..'

        bbox = draw.textbbox((0, 0), display_name, font=font)
        text_width = bbox[2] - bbox[0]
        x_pos = (size[0] - text_width) // 2
        draw.text((x_pos, y_pos), display_name, fill='#1a237e', font=font)

        y_pos += 16

        # Draw formula
        bbox = draw.textbbox((0, 0), formula, font=font_small)
        text_width = bbox[2] - bbox[0]
        x_pos = (size[0] - text_width) // 2
        draw.text((x_pos, y_pos), formula, fill='#424242', font=font_small)

        # Save
        combined.save(output_path, 'PNG')
        return True

    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def validate_structure(smiles: str, expected_formula: str) -> dict:
    """Validate a structure and return info."""
    result = {
        'valid': False,
        'rings': [],
        'carbons': 0,
        'mw': 0,
        'calculated_formula': '',
        'formula_match': False
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result

        result['valid'] = True
        result['mw'] = round(Descriptors.MolWt(mol), 2)
        result['carbons'] = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
        result['calculated_formula'] = rdMolDescriptors.CalcMolFormula(mol)

        ring_info = mol.GetRingInfo()
        result['rings'] = sorted([len(r) for r in ring_info.AtomRings()])

        # Check formula match
        result['formula_match'] = (
            result['calculated_formula'].replace(' ', '') ==
            expected_formula.replace(' ', '')
        )

    except Exception:
        pass

    return result


def main():
    print("=" * 80)
    print("GENERATING LABELED STRUCTURES FOR ALL 130 DROPDOWN COMPOUNDS")
    print("=" * 80)

    report = {
        'total': len(COMPOUNDS),
        'successful': 0,
        'failed': 0,
        'formula_mismatches': 0,
        'compounds': []
    }

    # Sort compounds alphabetically
    sorted_compounds = sorted(COMPOUNDS.keys())

    for i, name in enumerate(sorted_compounds, 1):
        compound = COMPOUNDS[name]
        print(f"\n[{i:3d}/{len(COMPOUNDS)}] {name}")

        # Validate
        validation = validate_structure(compound.smiles, compound.molecular_formula)

        print(f"  SMILES: {compound.smiles[:50]}{'...' if len(compound.smiles) > 50 else ''}")
        print(f"  Formula: {compound.molecular_formula} (calc: {validation['calculated_formula']})")
        print(f"  MW: {validation['mw']}")
        print(f"  Rings: {validation['rings']}")

        if not validation['formula_match']:
            print(f"  WARNING: Formula mismatch!")
            report['formula_mismatches'] += 1

        # Generate image
        output_path = os.path.join(OUTPUT_DIR, f"{name}.png")
        success = generate_labeled_structure(
            name,
            compound.smiles,
            compound.molecular_formula,
            output_path
        )

        if success:
            print(f"  Image: {output_path}")
            report['successful'] += 1
        else:
            print(f"  FAILED to generate image")
            report['failed'] += 1

        report['compounds'].append({
            'name': name,
            'smiles': compound.smiles,
            'formula_db': compound.molecular_formula,
            'formula_calc': validation['calculated_formula'],
            'mw': validation['mw'],
            'rings': validation['rings'],
            'valid': validation['valid'],
            'formula_match': validation['formula_match'],
            'image_generated': success
        })

    # Save report
    report_path = os.path.join(OUTPUT_DIR, 'validation_report.json')
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total compounds: {report['total']}")
    print(f"Images generated: {report['successful']}")
    print(f"Failed: {report['failed']}")
    print(f"Formula mismatches: {report['formula_mismatches']}")
    print(f"\nReport saved to: {report_path}")

    # Print all compounds
    print("\n" + "=" * 80)
    print("ALL 130 COMPOUNDS:")
    print("=" * 80)
    for i, c in enumerate(report['compounds'], 1):
        status = "✓" if c['image_generated'] else "✗"
        formula_status = "" if c['formula_match'] else " [FORMULA MISMATCH]"
        print(f"{i:3d}. {status} {c['name']:30s} {c['formula_db']:12s} rings={c['rings']}{formula_status}")

    return report


if __name__ == '__main__':
    main()
