import os
import json
import sys
import glob
sys.path.append(os.getcwd())
from gecko_web import reaction_tree

DATA_DIR = "/Users/rgpro/Desktop/GECKO/data/output"
JOB_STATE_FILE = "/Users/rgpro/Desktop/GECKO/data/job_state.json"

def regenerate_all():
    print("Regenerating reaction trees for all jobs...")
    count = 0
    if not os.path.exists(DATA_DIR):
        print("No output directory found.")
        return

    job_state = {}
    if os.path.exists(JOB_STATE_FILE):
        try:
            with open(JOB_STATE_FILE, "r") as f:
                job_state = json.load(f)
        except Exception as e:
            print(f"Error reading job_state.json: {e}")

    for job_id in os.listdir(DATA_DIR):
        job_dir = os.path.join(DATA_DIR, job_id)
        if not os.path.isdir(job_dir):
            continue
            
        print(f"Processing {job_id}...")
        try:
            voc_name = "unknown_voc"
            if job_id in job_state:
                voc_name = job_state[job_id].get("voc_name", "unknown_voc")
                print(f"  -> Found VOC: {voc_name}")

            # We pass the real voc_name so the fix works
            tree_data = reaction_tree.parse_reaction_tree(job_dir, voc_name)
            
            if tree_data and tree_data.get('nodes'):
                with open(os.path.join(job_dir, "reaction_tree.json"), "w") as f:
                    json.dump(tree_data, f)
                print(f"  -> Updated {len(tree_data['nodes'])} nodes.")
                count += 1
            else:
                print("  -> No nodes found or failed to parse.")
                
        except Exception as e:
            print(f"  -> Failed: {e}")

    print(f"Done. Updated {count} jobs.")

if __name__ == "__main__":
    regenerate_all()
