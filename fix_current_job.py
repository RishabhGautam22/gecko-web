import os
import sys
import json

# Add gecko_web to path
sys.path.append(os.path.join(os.getcwd(), "gecko_web"))
from gecko_web import reaction_tree

def fix_current_job():
    output_dir = "/Users/rgpro/Desktop/GECKO/data/output/boxmodel-e717ff7f49de" 
    voc_name = "isopentane"
    
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist.")
        return

    print(f"Regenerating tree for {output_dir}...")
    
    # Use the fixed reaction_tree.parse_reaction_tree
    tree_data = reaction_tree.parse_reaction_tree(output_dir, voc_name)
    
    # Check if fix worked
    fixed = False
    for node in tree_data.get('nodes', []):
        if node['id'] in ['C05000', 'G05000'] and node['smiles'] == 'CC(C)CC':
            fixed = True
            print(f"Verified node {node['id']} has correct SMILES: {node['smiles']}")
            
    if fixed:
        print("SUCCESS: Fix verified.")
        # Save to disk
        output_file = os.path.join(output_dir, "reaction_tree.json")
        with open(output_file, 'w') as f:
            json.dump(tree_data, f, indent=2)
        print(f"Saved updated tree to {output_file}")
    else:
        print("FAIL: Still not showing CC(C)CC")

if __name__ == "__main__":
    fix_current_job()
