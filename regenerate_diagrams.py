import os
import json
import sys
import glob

# Add current directory to path
sys.path.append(os.getcwd())
try:
    from gecko_web import pathway_visualizer
except ImportError:
    # If not running from root, try adding it
    sys.path.append(os.path.join(os.getcwd(), 'gecko_web'))
    from gecko_web import pathway_visualizer

DATA_DIR = "/Users/rgpro/Desktop/GECKO/data/output"
JOB_STATE_FILE = "/Users/rgpro/Desktop/GECKO/data/job_state.json"

def regenerate_diagrams():
    print("Regenerating diagrams for all jobs using fixed reaction trees...")
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
            
        json_path = os.path.join(job_dir, "reaction_tree.json")
        if not os.path.exists(json_path):
            print(f"Skipping {job_id} (no reaction_tree.json)")
            continue

        print(f"Processing {job_id}...")
        
        # Determine VOC name for title
        voc_name = job_id
        if job_id in job_state:
            voc_name = job_state[job_id].get("voc_name", job_id)
        
        try:
            with open(json_path, 'r') as f:
                tree_data = json.load(f)
            
            if not tree_data.get('nodes'):
                print(f"  -> Skipping (empty nodes)")
                continue
                
            # Regenerate PNG
            output_base = os.path.join(job_dir, "pathway_diagram")
            try:
                pathway_visualizer.generate_pathway_diagram(
                    tree_data=tree_data,
                    output_path=output_base,
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    format='png'
                )
                print(f"  -> Generated PNG")
            except Exception as e:
                print(f"  -> Failed PNG: {e}")

            # Regenerate SVG
            try:
                pathway_visualizer.generate_pathway_diagram(
                    tree_data=tree_data,
                    output_path=output_base,
                    title=f"Oxidation Mechanism: {voc_name.upper()}",
                    format='svg'
                )
                print(f"  -> Generated SVG")
            except Exception as e:
                print(f"  -> Failed SVG: {e}")
                
            count += 1
                
        except Exception as e:
            print(f"  -> Failed: {e}")

    print(f"Done. Regenerated diagrams for {count} jobs.")

if __name__ == "__main__":
    regenerate_diagrams()
