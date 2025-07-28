import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def read_pul_file(filepath):
    pul_category_sequences = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split by ':' to separate PUL id and categories string
            pul_id, categories_str = line.split(":", 1)
            # Split categories by the arrow symbol and strip spaces
            categories = [cat.strip() for cat in categories_str.split("→")]
            pul_category_sequences[pul_id.strip()] = categories
    return pul_category_sequences


def load_pul_organisms(filepath):
    pul_to_organism = {}
    with open(filepath, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            pul_id, organism = line.strip().split(',', 1)
            pul_to_organism[pul_id] = organism
    return pul_to_organism


def plot_pul_categories(pul_category_sequences, pul_to_organism=None):
    # Assign each category a unique color
    categories = sorted({cat for seq in pul_category_sequences.values() for cat in seq})
    # color_map = {cat: plt.cm.tab20(i % 20) for i, cat in enumerate(categories)}
    import seaborn as sns
    palette = sns.color_palette("Set2", n_colors=len(categories))
    color_map = {cat: palette[i] for i, cat in enumerate(categories)}


    # pul_labels = []
    # for pul_id in pul_category_sequences:
    #     org = pul_to_organism.get(pul_id, "") if pul_to_organism else ""
    #     label = f"{pul_id} ({org})" if org else pul_id
    #     pul_labels.append(label)


    fig, ax = plt.subplots(figsize=(12, len(pul_category_sequences) * 0.5))

    pul_ids = list(pul_category_sequences.keys())

    for i, (pul_id, seq) in enumerate(pul_category_sequences.items()):
        # Sort each PUL's gene categories according to the order in categories
        # seq = sorted(seq, key=lambda x: categories.index(x))
        for j, category in enumerate(seq):
            ax.barh(i, 1, left=j, color=color_map[category], edgecolor='black')
            
        if pul_to_organism and pul_id in pul_to_organism:
            ax.text(
                x=len(seq) + 0.5, y=i,
                va='center', ha='left',
                fontsize=8, color='gray',  s=r"$\mathit{" + pul_to_organism[pul_id].replace(' ', r'\ ') + "}$"
            )

    # Customize ticks
    ax.set_yticks(range(len(pul_category_sequences)))
    ax.set_yticklabels(pul_ids) # list(pul_category_sequences.keys())
    ax.set_xlabel("Gene Position")
    ax.set_title("PUL Category Sequences (Sorted by Length)")

    # Legend
    handles = [mpatches.Patch(color=color_map[cat], label=cat) for cat in categories]
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title="Category")

    plt.tight_layout()
    plt.show()
    
# def plot_pul_categories(pul_category_sequences, pul_to_organism=None):
#     import matplotlib.cm as cm
#     custom_order = ["CAZyme", "Transporter", "SusC", "SusD", "STP", "TF-STP", "Sulfatase", "Peptidase", "Hypothetical", "Other"] # "TF"

#     # Filter actual categories used in data (preserve order from custom_order)
#     used_categories = [cat for cat in custom_order if any(cat in seq for seq in pul_category_sequences.values())]

#     # Assign consistent colors
#     cmap = cm.get_cmap("Set3", len(used_categories))
#     color_map = {
#         "CAZyme": "#1f77b4",       # blue
#         "Transporter": "#ff7f0e",  # orange
#         "SusC": "#2ca02c",         # green
#         "SusD": "#2ca02c",         # green
#         "STP": "#d62728",          # red
#         # "TF": "#9467bd",           # purple
#         "TF-STP": "#9467bd",      
#         "Sulfatase": "#8c564b",    # brown
#         "Peptidase": "#e377c2",    # pink
#         "Hypothetical": "#17becf" , # cyan
#         "Other": "#7f7f7f"        # gray
#     }

#     # Sort PULs from largest to smallest by sequence length
#     sorted_puls = sorted(pul_category_sequences.items(), key=lambda x: len(x[1]), reverse=True)

#     # Plot
#     fig, ax = plt.subplots(figsize=(12, len(sorted_puls) * 0.5))

#     for i, (pul_id, seq) in enumerate(sorted_puls):
#         # Sort each PUL's gene categories according to custom_order
#         sorted_seq = sorted(seq, key=lambda x: custom_order.index(x) if x in custom_order else 999)
#         for j, category in enumerate(sorted_seq):
#             ax.barh(i, 1, left=j, color=color_map.get(category, "gray"), edgecolor='black')
            
#         if pul_to_organism and pul_id in pul_to_organism:
#             ax.text(
#                 x=len(seq) + 0.3, y=i,
#                 va='center', ha='left',
#                 fontsize=8, color='gray',  s=r"$\mathit{" + pul_to_organism[pul_id].replace(' ', r'\ ') + "}$")

#     # Customize ticks
#     ax.set_yticks(range(len(sorted_puls)))
#     ax.set_yticklabels([pul_id for pul_id, _ in sorted_puls])
#     ax.set_xlabel("Gene Position (Grouped by Custom Category Order)")
#     ax.set_title("PUL Category Sequences (Grouped by Custom Order, Largest to Smallest)")
#     ax.invert_yaxis()

#     # Legend
#     handles = [mpatches.Patch(color=color_map[cat], label=cat) for cat in used_categories]
#     plt.legend(handles=handles, title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')

#     plt.tight_layout()
#     plt.show()

if __name__ == "__main__":
    filename = "xylan_pul_category_sequences_with_susCD.txt"  # Change this to your filename
    pul_to_organism = load_pul_organisms("xylan_pul_org.csv")
    pul_category_sequences = read_pul_file(filename)
    plot_pul_categories(pul_category_sequences, pul_to_organism=pul_to_organism)
