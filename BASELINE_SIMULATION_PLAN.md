# Baseline Simulation & Synthetic Data Fallback Plan

## Overview
To ensure the application always presents valid, scientifically rigorous visualizations—even when user-initiated simulations fail or return partial data—we will implement a **Baseline Simulation Fallback System**. 

This system ensures "no placeholders" by substituting missing real data with explicitly labeled "Baseline Data" (e.g., from a verified Alpha-Pinene reference run) or "Synthetic Data" (statistically derived).

## 1. Hierarchy of Data Sources
The system will attempt to load data in this order of preference:

1.  **Real Simulation Data:** Actual `concentrations.nc` or `.out` files from the current Job ID.
2.  **Cached Reference Data:** Pre-computed results for the specific VOC if available in the library.
3.  **Baseline Surrogate Data:** A verified dataset from a similar class of compound (e.g., Alpha-Pinene for terpenes, Toluene for aromatics), scaled to match the requested initial conditions.
4.  **Synthetic Fallback:** Statistically generated time series based on mechanism properties (already implemented in `postprocessing.py`).

## 2. Implementation Logic

### A. Detection Phase (in `combined_workflow.py` & `postprocessing.py`)
- **Action:** Check for existence and validity (non-zero size, readable headers) of output files.
- **Trigger:** If `concentrations.nc` is missing or empty.
- **Response:**
    - Log the missing data event.
    - Set `job_state.data_source = "synthetic"` or `"baseline_surrogate"`.
    - Select the appropriate fallback dataset.

### B. Baseline Injection (New Module: `baseline_manager.py`)
- **Structure:**
  ```python
  def get_baseline_data(voc_class: str, initial_ppb: float):
      # Load reference dataset (e.g., verified APINENE runs)
      df = load_reference("apinene_high_nox")
      # Scale concentrations linearly by initial_ppb ratio
      scale_factor = initial_ppb / ref_initial_ppb
      return df * scale_factor
  ```
- **Metadata Injection:**
  - Add a flag to the DataFrame: `df.attrs['is_synthetic'] = True`.
  - Add a "Note" column describing the source.

### C. Visual Indication (Frontend)
- **UI Banner:** "⚠️ DISPLAYING BASELINE SURROGATE DATA" (Yellow/Orange).
- **Plot Overlay:** Watermark or subtitle on every plot: "Simulated Data (Baseline: Alpha-Pinene)".
- **Report Appendix:** The "Data Quality Appendix" will explicitly state:
  > "Primary simulation output was unavailable. Results presented are derived from [SOURCE] for demonstration purposes."

## 3. Automation Steps for "Auto-Run"
1.  **Pre-computation:** Run trusted reference simulations for 5 key classes (Terpene, Isoprene, Aromatic, Alkane, Alkene) at standard conditions (High/Low NOx).
2.  **Storage:** Store these results in `data/library/baselines/`.
3.  **Fallback Trigger:** Update `combined_workflow.py`:
    ```python
    if not job_success:
       logger.warn("Job failed, loading baseline...")
       baseline = select_baseline(voc_structure)
       copy_baseline_to_output(baseline, output_dir)
    ```

## 4. User Interaction (Pop-up Feature)
- **Problem:** Ambiguity when input parameters are missing.
- **Solution:** "Smart Wizard" modal.
    - If a user inputs a custom SMILES but no vapor pressure model:
    - **Popup:** "Vapor Pressure Model Required. Recommended: Nannoolal. Use default?"
    - **Action:** If user accepts, inject default params into `job_state.json`.

## 5. Next Steps
1.  Verify `postprocessing.py` correctly handles the `allow_synthetic` flag (Done).
2.  Add `is_synthetic` flag to `job_state.json` schema.
3.  Update frontend `app.js` to read this flag and show the banner (Partially done).
4.  Curate the `data/library/baselines/` dataset.
