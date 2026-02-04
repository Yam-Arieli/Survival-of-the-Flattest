import pandas as pd

# 1. Load your merged table
df = pd.read_csv("ABF1_with_Metadata.csv")

# 2. Check which column name actually exists
# In Table S1 it is 'Ecological origins', in some scripts we renamed it 'Origin_Type'
if 'Origin_Type' in df.columns:
    col_name = 'Origin_Type'
elif 'Ecological origins' in df.columns:
    col_name = 'Ecological origins'
else:
    print(f"Error: Could not find ecological column. Columns are: {df.columns.tolist()}")
    exit()

# 3. Define your ecological targets
# Note: 'Industrial' might be listed as 'Industrial/Domesticated' or similar in some versions
target_origins = ['Soil', 'Lab strain', 'Industrial']

# 4. Filter the data
# We filter for the selected origins AND ensure the Environment is YPD (Glucose-rich)
filtered_df = df[
    (df[col_name].isin(target_origins)) &
    (df['Environment'] == 'YPD')
]

# 5. Save the filtered result
filtered_df.to_csv("Filtered_ABF1_Glucose.csv", index=False)

# 6. Summary of the results
print(f"Total rows in subset: {len(filtered_df)}")
print("-" * 30)
print(f"Breakdown of filtered strains (using column '{col_name}'):")
print(filtered_df[col_name].value_counts())
