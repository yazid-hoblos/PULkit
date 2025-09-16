# Versatile Toolkit for PUL Systems Detection

This repository offers a wide range of scripts developed within the scope of a project aimed at **PUL** (Polysaccharide Utilization Loci) systems detection at the pangenome level. It encompasses three main broad themes: 

<ol type="1" style="list-style: none;">
  <li> PUL models formulation</li>
  <li> Helper scripts for PUL systems detection and validation</li>
  <li> Comparartive analysis of predictions across approaches</li>
</ol>

## 1. PUL Models Formulation **`/pul_models/`** 

PUL prediction models remain currently limited to two main approaches in the literature:

1. The _SusCD_ model introduced by _Terrapon et al._ [[1]](#ref1) in 2015, formally defining PULs as operon-like structures centered around a SusC-SusD adjacent pair (TonB-dependent transporter). It primarily concerns _Bacteroidetes_ species as it was defined and validated on them.

2. To further generalize PUL systems definition, dbCAN2 [[2]](#ref2) introduced _CAZyme Gene Clusters (CGCs)_ in 2018 as clusters of at least one CAZyme, one transporter component (TC), and one transcription factor (TF), separated by a maximum of two non-signature genes.

Here, the potential of formulating specific PUL models in a more rigorous manner based on experimentally-validated PULs is explored, primarily on the basis of substrate-based classification. The main challenge remains the limited availability of curated PUL data. One potential way to circumvent this issue is making use of the predictions obtained based on more general models, such as the two highlighted. 

To this end, high-confidence predictions across diverse species and genera, preferably at the pangenome level, could be analyzed to reveal distinctive signature patterns. Adopting this approach within the context of a specific inquiry, e.g. PULs associated with a specific substrate (which could be predicted for unvalidated PULs), would increase the chances of uncovering insightful patterns to guide new models formulation. Thus, the suggested approach is bidirectional: current models are used to generate predictions, which guide the iterative refinement of the models. 

### 1.1. Exploring PUL Signature Patterns   
**`/pul_models/signature_patterns/`**

An initial attempt to explore patterns across heterogenous experimentally-validated PULs could be found in this directory, primarily focused on: 
- Component pattern analysis 
- Signature motifs within PULs
- Statistical patterns across PULs

A detailed description of the analysis that was done could be found in the directory's corresponding report: [README](pul_models/signature_patterns/README.md).

### 1.2. Curated PUL Data 
**`/pul_models/pul_data/`**

- Curated empirically-validated PUL systems
- _SusCD_ HMM profiles
- The scripts used for data accession and analysis are in [bin/](pul_models/bin/)

More details on the scripts are covered step-by-step in the README of **substrate-specific/**: [README](pul_models/substrate-specific/README.md).

### 1.3. Substrate-specific PUL Models 
**`/pul_models/substrate-specific/`**

- Models formultion for PULs associated to one substrate
- Curated data and initial findings for xylan PULs

More details in associated report: [README](pul_models/substrate-specific/README.md).

## 2. Utils for PUL Systems Detection & Validation


### 📁 `/annotation/` Pangenome signature gene families annotation 

- `write_metadata_to_pangenomes.sh` - Add metadata to pangenome files
- `set_hmm_thresholds.sh` - Configure HMM profile thresholds
- `extract_protein_families_per_pangenome.sh` - Extract protein family assignments
- `extract_metadata_per_pangenome.sh` - Extract pangenome metadata
- `diamond_tp_annotation_per_pangenome.sh` - DIAMOND-based transporter annotation
- `dbcan_cazymes_annotation.sh` - dbCAN CAZyme annotation pipeline
- `combine_cazyme_metadata.sh` - Merge CAZyme annotation results
- `add_tp_metadata.sh` - Add transporter metadata to datasets
- `add_tf_stp_metadata.sh` - Add transcription factor and STP metadata

### 📁 `/panorama_utils/` PANORAMA predictions analysis

- `systems_overlap_heatmap.py` - Generate overlap heatmaps between system predictions
- `reduce_projected_system_ids.py` - Simplify system ID mappings
- `process_systems_overlaps.py` - Process overlap statistics between systems
- `find_unmatched_systems.py` - Identify systems without matches
- `find_exact_system_matches.py` - Find identical system architectures
- `compare_systems_overlap.py` - Calculate overlap metrics between system sets
- `compare_panorama_systems.py` - Compare Panorama predictions across datasets

### 📁 `/extract_pul_data/` PUL data extraction & manipulation for predictions validation

- `map_cazy_old_tags.py` - Map legacy CAZy identifiers to current format
- `get_curated_cazy_puls.py` - Extract curated PULs from CAZy database
- `fetch_taxon_strains.py` - Retrieve strain information by taxonomy
- `fetch_cazy_puls_per_species.py` - Species-specific PUL extraction
- `examine_small_dbcan_puls.py` - Analyze small/minimal PULs from dbCAN
- `draw_pul_components_frequency_dist.py` - Visualize component frequency distributions
- `cazy_puls_size_dist.py` - Analyze PUL size distributions from CAZy

### 📁 `/utils/` Diverse helper scripts

- `get_genomes_per_species.sh` - Extract genome datasets for specific species
- `extract_ids_mapping.sh` - Create ID mapping files between databases
- `create_gff_tsv.sh` - Convert GFF files to TSV format for analysis
- `create_fasta_tsv.sh` - Process FASTA files into tabular format
- `compare_panorama_cgc_finder_predictions.sh` - Compare CGC predictions between tools
- `cgc_finder_across_pangenomes.sh` - Run CGC finder on multiple pangenomes
- `strains_per_species_table.sh` - Generate strain count tables by species
- `cazy_puls_curation.sh` - Curate PUL data from CAZy database
- `archived_cgc_finder_across_pangenome.sh` - Legacy CGC analysis script


## 3. Comparative Analysis of PANORAMA, dbCAN, and CAZy Predictions 

### **`/comparative_analysis/`**

- `umap_dbscan_systems.py` - UMAP dimensionality reduction with DBSCAN clustering
- `systems_similarity_heatmap.py` - Create similarity heatmaps between systems
- `venn_diagram.py` - Generate Venn diagrams for dataset overlaps
- `systems_pca_clustering.py` - Principal component analysis and clustering
- `systems_overlap_scatter_plots.py` - Scatter plot visualizations of system overlaps
- `systems_overlap_coverage_plot.py` - Coverage analysis plots
- `plot_pul_partitions.py` - Visualize PUL partitioning results
- `panorama_cazy_barplot.py` - Bar plots comparing Panorama vs CAZy predictions
- `pan_genus_umap.py` - Genus-level UMAP projections
- `pan_genus_heatmaps.py` - Genus-level comparative heatmaps
- `match_panorama_cazy_dbcan_systems.py` - Match systems across prediction tools
- `completeness_plot.py` - Visualize dataset completeness metrics
- `combine_images.py` - Combine multiple plots into publication figures
- `cgc_finder_predictions_dist.py` - Distribution analysis of CGC predictions

### `/comparative_analysis/outputs/` 

- Generated visualizations and figures
- Some intermediate processing files


## Typical Workflow:

1. **PUL Models Formulation** (`/pul_models/`)
   
   ↓ Explore signature patterns and generate models

2. **Annotation** (`/annotation/`)
   
   ↓ Annotate signature gene families in pangenomes

3. **PUL Detection**
   
   ↓ Run PANORAMA to detect PULs in pangenomes

4. **Pangenome-level Insights** (`/panorama_utils/`)
   
   ↓ Examine PANORAMA predictions

5. **Utilities** (`/utils/`)
   
   ↓ Support scripts for data processing and format conversion

6. **Predicted and Empirical PULs Extraction** (`/extract_pul_data/`)
   
   ↓ Validation of predictions with empirical PULs and other approaches predictions

7. **Comparative Analysis** (`/comparative_analysis/`)
   
   ↓ Compare across species and approaches


## References

1. <a id="ref1"></a>Nicolas Terrapon, Vincent Lombard, Harry J. Gilbert, and Bernard Henrissat. Automatic prediction of polysaccharide utilization loci in bacteroidetes species. Bioinformatics, 31(5):647–655, March 2015. doi: [10.1093/bioinformatics/btu716](https://doi.org/10.1093/bioinformatics/btu716).

2. <a id="ref2"></a>Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K. Busk, Ying Xu, and Yanbin Yin. dbcan2: a meta server for automated carbohydrate-active enzyme annotation. Nucleic Acids Research, 46(W1):W95–W101, 2018. doi: [10.1093/nar/gky418](https://doi.org/10.1093/nar/gky418).