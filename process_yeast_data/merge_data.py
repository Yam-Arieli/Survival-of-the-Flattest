import pandas as pd
import io

# 1. Load your sequences
sequences_df = pd.read_csv("WT_ABF1_Sequences.csv")

def load_clean_metadata(filename):
    # Read the file as raw text first to find where the real table starts
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find the line that actually contains our headers
    start_row = 0
    for i, line in enumerate(lines):
        if "Standardized name" in line:
            start_row = i
            break
    
    # Reload the CSV starting from that specific row
    df = pd.read_csv(filename, skiprows=start_row)
    # Clean up column names (removes spaces/newlines)
    df.columns = df.columns.str.strip()
    return df

# 2. Load the metadata CSV using our cleaner function
metadata_df = load_clean_metadata("TableS1_metadata.csv")

# 3. Check and Merge
required_cols = ['Standardized name', 'Ecological origins', 'Geographical origins']
missing = [c for c in required_cols if c not in metadata_df.columns]

if missing:
    print(f"Error: Still missing columns: {missing}")
    print(f"Pandas found these columns instead: {metadata_df.columns.tolist()}")
else:
    # Merge on 'Strain' and 'Standardized name'
    merged_df = pd.merge(
        sequences_df,
        metadata_df[required_cols],
        left_on='Strain',
        right_on='Standardized name',
        how='left'
    )
    
    # Drop duplicate and save
    merged_df.drop(columns=['Standardized name'], inplace=True)
    merged_df.to_csv("ABF1_with_Metadata.csv", index=False)
    print("Success! Created 'CDC36_with_Metadata.csv' with ecological data.")
