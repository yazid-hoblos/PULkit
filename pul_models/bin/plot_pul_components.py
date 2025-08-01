'''Plots a bar plot representation of PUL category sequences.'''

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


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


def plot_pul_categories(pul_category_sequences, output_file, pul_to_organism_file=None,
                        custom_order=None, sort_by_length=True):

    if pul_to_organism_file:
        pul_to_organism = load_pul_organisms(pul_to_organism_file)

    # Sort PULs by length if specified
    pul_items = sorted(pul_category_sequences.items(), key=lambda x: len(x[1])) if sort_by_length else pul_category_sequences.items()

    # Get categories
    if custom_order:
        used_categories = [cat for cat in custom_order if any(cat in seq for _, seq in pul_items)]
        # convert all categories not in custom_order to "Other"
        pul_items = [(pul_id, [cat if cat in used_categories else "Other" for cat in seq]) for pul_id, seq in pul_items]
        color_palette = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#17becf", "#7f7f7f"
        ]
        color_map = {cat: color_palette[i % len(color_palette)] for i, cat in enumerate(used_categories)}
    else:
        used_categories = sorted({cat for _, seq in pul_items for cat in seq})
        palette = sns.color_palette("Set2", n_colors=len(used_categories))
        color_map = {cat: palette[i] for i, cat in enumerate(used_categories)}

    # Plot setup
    fig, ax = plt.subplots(figsize=(20, len(pul_items) * 0.35))
    
    for i, (pul_id, seq) in enumerate(pul_items):
        categories_seq = (
            sorted(seq, key=lambda x: custom_order.index(x) if x in custom_order else 999)
            if custom_order else seq
        )

        for j, category in enumerate(categories_seq):
            ax.barh(i, 1, left=j, color=color_map.get(category, "gray"), edgecolor='black')

        if pul_to_organism and pul_id in pul_to_organism:
            ax.text(x=len(seq) + 0.3, y=i, va='center', ha='left',
                    fontsize=10, color='gray',
                    s=r"$\mathit{" + pul_to_organism[pul_id].replace(' ', r'\ ') + "}$")

    # Ticks and labels
    ax.set_yticks(range(len(pul_items)))
    ax.set_yticklabels([pul_id for pul_id, _ in pul_items])
    ax.set_xlabel("Gene Position")
    ax.set_title("PUL Category Sequences" + 
                 (" (Grouped by Custom Order)" if custom_order else ""))
    ax.invert_yaxis()

    # Legend
    handles = [mpatches.Patch(color=color_map[cat], label=cat) for cat in used_categories]
    plt.legend(handles=handles, title="Category", bbox_to_anchor=(1, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f"pul_component_patterns/{output_file}", dpi=500)



if __name__ == "__main__":
    filename = "pul_component_patterns/category_sequences.txt"  
    pul_category_sequences = read_pul_file(filename)
    org_file = 'data/xylan_pul_org.csv'
    plot_pul_categories(pul_category_sequences, 'all_categories_considered.png', pul_to_organism_file=org_file)
    custom_order = ["CAZyme", "Transporter", "SusC", "SusD", "STP", "TF-STP", "Sulfatase", "Peptidase", "Hypothetical", "Other"] # all TFs were found to be STPs as well
    plot_pul_categories(pul_category_sequences, 'custom_ordered_plot.png', pul_to_organism_file=org_file, custom_order=custom_order)

