import os
import pandas as pd
import xarray as xr
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def inspect_data():
    output_dir = "/Users/rgpro/Documents/GEKO/data/output/test-job-2"
    nc_file = os.path.join(output_dir, "alpha-pinene.nc")
    dict_file = os.path.join(output_dir, "dictionary.out")
    
    print(f"Checking files in {output_dir}")
    
    # 1. Inspect NetCDF keys
    if os.path.exists(nc_file):
        try:
            ds = xr.open_dataset(nc_file)
            df = ds.to_dataframe().reset_index()
            ds.close()
            
            print("\n--- NetCDF Columns (First 20) ---")
            print(df.columns[:20].tolist())
            
            if 'Species' in df.columns:
                print("\n--- Species Column Sample ---")
                sample = df['Species'].head(10).tolist()
                print(sample)
                
                # Decode and strip
                decoded = [s.decode('utf-8').strip() for s in sample]
                print(f"Decoded: {decoded}")
                
                # Check for 'G' prefix
                stripped = [s[1:] if s.startswith('G') else s for s in decoded]
                print(f"Stripped 'G': {stripped}")
                
                # Check matches
                if os.path.exists(dict_file):
                    from gecko_web.postprocessing import parse_dictionary
                    meta = parse_dictionary(dict_file)
                    dict_ids = set(meta['id'])
                    
                    matches = [s for s in stripped if s in dict_ids]
                    print(f"Matches found in sample: {len(matches)} / {len(sample)}")
                    print(f"Matching IDs: {matches}")
            
            # Check if 'Species' column exists or if species are columns
            # Usually in these outputs, species are columns.
            # But postprocessing.py expects a 'Species' column in a long-format dataframe?
            # Let's check the shape.
            print(f"\nDataFrame Shape: {df.shape}")
            
            # If species are columns, we need to melt it?
            # postprocessing.py:
            # merged = pd.merge(data_df, species_meta, left_on='Species', right_on='code', how='left')
            # This implies data_df has a 'Species' column.
            
            if 'Species' not in df.columns:
                print("\n'Species' column NOT found. Species are likely columns.")
                # Sample column names that might be species
                potential_species = [c for c in df.columns if c not in ['time', 'Time', 'id', 'index']]
                print(f"Sample potential species columns: {potential_species[:10]}")
        except Exception as e:
            print(f"Error reading NetCDF: {e}")
    else:
        print("NetCDF file not found.")

    # 2. Inspect Dictionary keys
    if os.path.exists(dict_file):
        try:
            from gecko_web.postprocessing import parse_dictionary
            meta = parse_dictionary(dict_file)
            print("\n--- Dictionary Entries (First 10) ---")
            if not meta.empty:
                print(meta[['id', 'code', 'formula_str']].head(10))
            else:
                print("Dictionary parsed as empty.")
        except Exception as e:
            print(f"Error reading Dictionary: {e}")
    else:
        print("Dictionary file not found.")

if __name__ == "__main__":
    inspect_data()
