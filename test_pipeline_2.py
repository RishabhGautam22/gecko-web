import os
import sys
import logging
import time
from pydantic import BaseModel
from typing import Optional

# Mock classes from main.py
class ScenarioParams(BaseModel):
    nox_regime: str = "high"
    temperature_k: float = 298.0
    rh_percent: float = 50.0
    seed_aerosol_ug_m3: float = 10.0
    dilution_rate_s1: float = 0.0

class JobRequest(BaseModel):
    voc_name: str
    job_type: str
    scenario: Optional[ScenarioParams] = None

# Set environment variables
cwd = os.getcwd()
os.environ["GECKO_SOURCE_DIR"] = os.path.join(cwd, "docker/gecko_source")
os.environ["BOXMODEL_SOURCE_DIR"] = os.path.join(cwd, "docker/boxmodel_source")
os.environ["DATA_DIR"] = os.path.join(cwd, "data")

# Import the logic from gecko_web.main
# We need to add the current directory to sys.path to import gecko_web
sys.path.append(cwd)

try:
    from gecko_web.main import run_process, jobs
except ImportError as e:
    print(f"Error importing gecko_web.main: {e}")
    sys.exit(1)

# Configure logging to see output
logging.basicConfig(level=logging.INFO)

def test_pipeline():
    job_id = "test-job-2"
    voc_name = "alpha-pinene" # Use a known VOC from the mapping
    job_req = JobRequest(voc_name=voc_name, job_type="boxmodel")
    
    # Initialize job entry in the global jobs dict (as run_process expects it)
    jobs[job_id] = {
        "id": job_id,
        "voc_name": voc_name,
        "type": "boxmodel",
        "status": "queued",
        "logs": [],
        "created_at": time.time()
    }
    
    print(f"Starting test job {job_id} for {voc_name}...")
    try:
        run_process(job_id, job_req)
    except Exception as e:
        print(f"run_process raised exception: {e}")
    
    job = jobs[job_id]
    print(f"Job Status: {job['status']}")
    if job['status'] == 'failed':
        print(f"Error: {job.get('error')}")
        print("Logs:")
        for log in job['logs']:
            print(log)
    else:
        print("Success!")
        print("Logs:")
        for log in job['logs']:
            print(log)
        
        # Check if plots exist
        output_dir = os.path.join(os.environ["DATA_DIR"], "output", job_id)
        import glob
        plots = glob.glob(os.path.join(output_dir, "*.png"))
        print(f"Plots found: {len(plots)}")
        for p in plots:
            print(f" - {os.path.basename(p)}")

if __name__ == "__main__":
    test_pipeline()
