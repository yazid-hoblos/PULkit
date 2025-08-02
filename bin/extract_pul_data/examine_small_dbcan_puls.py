import pandas as pd
import sys 


if len(sys.argv) != 2:
    print("Usage: python draw_PULs.py <dbcan_pul_components.csv> [output_file.csv]")
    sys.exit(1)

puls_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else "puls_with_3_or_less_components.csv"

df = pd.read_csv(puls_file, sep=',')

# Count number of components per PUL
component_counts = df["PUL_ID"].value_counts()

# Find PULs with 4 or fewer components
small_puls = component_counts[component_counts <= 3].index

# Filter original dataframe for those PULs
filtered_df = df[df["PUL_ID"].isin(small_puls)]

# Save or inspect
filtered_df.to_csv(output_file, index=False)