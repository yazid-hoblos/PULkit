import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load data
df = pd.read_csv("comparisons/previous_matches.tsv", sep="\t", header=None, names=["Genome", "PUL_genes", "Mapped_proteins", "Matches_other", "Other_tool_len"])

# Compute percentages
df["%mapped"] = df["Mapped_proteins"] / df["PUL_genes"] * 100
df["%matched_in_tool"] = df["Matches_other"] / df["PUL_genes"] * 100
df["%coverage_of_tool"] = df["Matches_other"] / df["Other_tool_len"] * 100

# Set plot style
sns.set(style="whitegrid")

# Melt dataframe for easier plotting
df_melted = df.melt(id_vars=["Genome"], value_vars=["%mapped", "%matched_in_tool", "%coverage_of_tool"],
                    var_name="Metric", value_name="Percentage")

# Create line frequency plot
plt.figure(figsize=(10, 6))
sns.histplot(data=df_melted, x="Percentage", hue="Metric", element="step", stat="density", common_norm=False, bins=50)
plt.title("Frequency Distribution of Annotation Metrics")
plt.xlabel("Percentage")
plt.ylabel("Density")
plt.xlim(0, 100)
plt.tight_layout()
plt.savefig("outputs/previous_percentage_lineplot.png")
plt.close()
