import pandas as pd
from pathlib import Path

# --- CONFIGURATION ---
SOURCE_DIR = Path('./Data_APR_Validation')
OUTPUT_DIR = Path('./out_parquets')
SHEET_INDEX = 3

# Map Excel column names to the names used in validation loader script
# k = Excel Column, v = Script Column
COLUMN_MAP = {
    'Relative Time(h:min:s.ms)': 'time',
    'Voltage(V)': 'voltage',
    'Current(A)': 'current',
    'Step_Status': 'status'
}

def excel2parquet():
    # Create output directory if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find all xlsx files
    files = list(SOURCE_DIR.glob('*.xlsx'))
    
    if not files:
        print(f"No files found in {SOURCE_DIR}")
        return

    for file_path in files:
        print(f"Processing: {file_path.name}...")
        
        try:
            df = pd.read_excel(file_path, sheet_name=SHEET_INDEX)
            df = df.rename(columns=COLUMN_MAP)

            df.time = pd.to_timedelta(df.time).dt.total_seconds()

            # fix the problem where timer is periodically reset to 0
            diffs = df['time'].diff().fillna(0)
            diffs[diffs < 0] = 0 
            df['time'] = diffs.cumsum()

            # Ensure essential columns exist before saving
            required_cols = ['time', 'voltage', 'current', 'status']
            missing = [c for c in required_cols if c not in df.columns]
            
            if missing:
                print(f"Skipping {file_path.name}: Missing columns {missing}")
                continue

            # 4. Save to Parquet
            output_filename = file_path.stem + ".parquet"
            df.to_parquet(OUTPUT_DIR / output_filename, index=False)
            print(f"Successfully converted to {output_filename}")

        except Exception as e:
            print(f"Error processing {file_path.name}: {e}")

if __name__ == "__main__":
    excel2parquet()