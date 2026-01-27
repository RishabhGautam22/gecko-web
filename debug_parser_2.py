
import sys
import os
sys.path.append(os.getcwd())

from gecko_web.reaction_tree import gecko_to_smiles, KNOWN_SPECIES
from rdkit import Chem

formulas_to_test = [
    ("2N5005", "CH3CH(ONO2)CH(CH3)CH2(OO.)"),
    ("AH05000", ""), # Should test the fallback if we could mock dictionary, but direct formula is empty
    ("H05000", "CH3CH(OOH)CH(CH3)CH3"),
    ("Isopentane_Parent", "CH3CH2CH(CH3)CH3"),  # C05000 / G05000
    ("Generic_C5", "C5H12")
]

print("--- Testing Formula Parsing ---")
for code, formula in formulas_to_test:
    print(f"\nTesting {code}: '{formula}'")
    smiles = gecko_to_smiles(formula, code)
    print(f"  -> SMILES: '{smiles}'")
    
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print("  -> RDKit Valid: Yes")
        else:
            print("  -> RDKit Valid: NO")
    else:
        print("  -> Result Empty")

print("\n--- Checking specific dictionary entries ---")
# Manually check if I can simulate the 'Isopentane' override condition
# This requires me to verify the logic in gecko_to_smiles or mechanism_diagram
