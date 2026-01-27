import pandas as pd
import io

def parse_dictionary_old(lines):
    species_data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 14:
            continue
        try:
            entry = {
                'id': parts[0],
                'mw': float(parts[1]),
                'nC': int(parts[3]),
                'nH': int(parts[4]),
                'nN': int(parts[5]),
                'nO': int(parts[6]),
                'nS': int(parts[7]),
                'code': parts[12],
                'formula_str': parts[13]
            }
            species_data.append(entry)
        except ValueError:
            continue
    return pd.DataFrame(species_data)

def parse_dictionary_new(lines):
    species_data = []
    for line in lines:
        parts = line.strip().split()
        # We expect at least: ID ... MW ir C H N O S F Br Cl
        # That is 1 + ... + 1 + 9 = 11 parts minimum.
        if len(parts) < 11:
            continue
        
        try:
            # The last 9 are atoms
            # The one before that is MW
            # Everything before that is the ID/Name info
            
            # Atoms: ir, C, H, N, O, S, F, Br, Cl
            # Indices from end:
            # -1: Cl
            # -2: Br
            # -3: F
            # -4: S
            # -5: O
            # -6: N
            # -7: H
            # -8: C
            # -9: ir
            # -10: MW
            
            mw = float(parts[-10])
            ir = int(parts[-9])
            nC = int(parts[-8])
            nH = int(parts[-7])
            nN = int(parts[-6])
            nO = int(parts[-5])
            nS = int(parts[-4])
            
            # ID is the first part
            # Code/Formula might be in the middle
            # Based on log: "TE4001   C1H2-O-C1(CH3)CHO   TDE"
            # parts[0] = TE4001
            # parts[1] = C1H2-O-C1(CH3)CHO
            # parts[2] = TDE
            
            entry = {
                'id': parts[0],
                'mw': mw,
                'nC': nC,
                'nH': nH,
                'nN': nN,
                'nO': nO,
                'nS': nS,
                'code': parts[0], # Use ID as code for matching if needed
                'formula_str': parts[1] if len(parts) > 10 else ""
            }
            species_data.append(entry)
        except ValueError:
            continue
    return pd.DataFrame(species_data)

log_line = "TE4001   C1H2-O-C1(CH3)CHO                                                                                                         TDE             86.1  0  4  6  0  2  0  0  0  0"
lines = [log_line]

print("--- Old Parser ---")
df_old = parse_dictionary_old(lines)
print(df_old)

print("\n--- New Parser ---")
df_new = parse_dictionary_new(lines)
print(df_new)
