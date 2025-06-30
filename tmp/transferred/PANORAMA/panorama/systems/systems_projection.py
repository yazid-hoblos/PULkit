#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to project systems onto genomes.
"""

# default libraries
from __future__ import annotations

import itertools
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import logging
from typing import Dict, List, Set, Tuple
from multiprocessing import Lock, get_context
from pathlib import Path
import time

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism, Gene

# local libraries
from panorama.utils import mkdir, init_lock, conciliate_partition
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.utils import (get_metadata_to_families, dict_families_context,
                                    get_gfs_matrix_combination, check_needed_families)
from panorama.systems.system import System, SystemUnit
from panorama.systems.models import Family


# def has_short_path(graph: nx.Graph, node_list, n):
#     """
#     Checks if there exists at least one path of length less than `n`
#     connecting any two nodes in the given list of nodes in the graph.

#     Args:
#         graph (nx.Graph): the graph to search paths
#         node_list (list): List of nodes to check for paths.
#         n (int): The maximum length of the path to consider.

#     Returns:
#         bool: True if there exists at least one path of length less than `n`
#               connecting any two nodes in the list, False otherwise.
#     """
#     path_length = defaultdict(dict)
#     has_path = {node: False for node in node_list}
#     for i, node1 in enumerate(node_list):
#         for node2 in node_list[i + 1:]:
#             if not has_path[node2]:
#                 try:
#                     path_length[node1][node2] = nx.shortest_path_length(graph, source=node1, target=node2)
#                     if path_length[node1][node2] <= n:
#                         has_path[node1] = True
#                         has_path[node2] = True
#                         break
#                 except nx.NetworkXNoPath:
#                     continue
#     return all(has_path.values())

# def has_short_path(graph: nx.Graph, node_list, n):
#     """
#     Checks if all nodes in the list can be connected through paths of length <= n.
    
#     Returns:
#         bool: True if all nodes have short paths to at least one other node, False otherwise.
#     """
#     if len(node_list) <= 1:
#         return True
        
#     has_path = {node: False for node in node_list}
    
#     for i, node1 in enumerate(node_list):
#         for node2 in node_list[i + 1:]:
#             try:
#                 path_length = nx.shortest_path_length(graph, source=node1, target=node2)
#                 if path_length <= n:
#                     has_path[node1] = True
#                     has_path[node2] = True
#             except nx.NetworkXNoPath:
#                 continue
    
#     return all(has_path.values())

    
def has_short_path(graph: nx.Graph, node_list, n):
    """
    Checks whether all nodes in `node_list` are connected to each other
    with a path of length ≤ `n`.

    Args:
        graph (nx.Graph): the graph containing the nodes
        node_list (list): list of nodes to check
        n (int): max allowed path length between any pair

    Returns:
        bool: True if all node pairs are at most `n` apart, False otherwise
    """
    # sort node_list by their gene.position
    node_list = sorted(node_list, key=lambda x: x.position)
    # print(f"Checking nodes: {[n.position for n in node_list]} with n={n}")

    for i, u in enumerate(node_list[:-1]):
        if nx.shortest_path_length(graph, source=u, target=node_list[i+1]) > n:
            return False

    # for i, u in enumerate(node_list):
        # for v in node_list[i+1:]:
            # try:
            # if nx.shortest_path_length(graph, source=u, target=v) > n:
                # return False
            # except nx.NetworkXNoPath:
            #     return False
    return True


def project_unit_on_organisms(graph: nx.Graph, unit: SystemUnit, organism: Organism,
                              model_genes: Set[Gene], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                              association: List[str] = None) -> Tuple[List[List[str]], List[int], str]:
    """
    Project a system onto an organism's pangenome.

    Args:
        graph (nx.Graph): Genomic context graph of the system for the given organism.
        unit (SystemUnit): The unit to be projected.
        organism (Organism): The organism on which the system is to be projected.
        gene_fam2mod_fam (Dict[str, Set[Family]]): A dictionary mapping gene families to model families.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[List[List[str]], List[int], str]: A tuple containing:
            - A list of projected system information for the organism.
            - A list with counts of each system organization type (strict, extended, split).
            - The reconciled system partition.
    """

    def write_projection_line(gene: Gene) -> List[str]:
        """
        Write the projection information for a single gene.

        Args:
            gene (Gene): Gene to write projection for.

        Returns:
            List[str]: List of elements representing the projection for the gene.
        """
        partition = gene.family.named_partition 
        line_projection = [gene.family.name, partition, fam_annot, fam_sec, gene.ID,
                           gene.local_identifier, gene.contig.name, gene.start, gene.stop, gene.strand, gene.is_fragment,
                           sys_state_in_org, gene.product]

        # if 'RGPs' in association:
        #     rgp = gene.RGP
        #     if rgp is not None:
        #         unit.add_region(rgp)
        #         line_projection.append(str(rgp))
        #     else:
        #         line_projection.append('')

        # if 'spots' in association:
        #     spot = gene.spot
        #     if spot is not None:
        #         unit.add_spot(gene.spot)
        #         line_projection.append(str(spot))
        #     else:
        #         line_projection.append('')

        return list(map(str, [unit.name, sub_id, organism.name] + line_projection))

    projection = []
    partitions = set()
    counter = [0, 0, 0]  # count strict, extended, and split CC

    sub_id = 1
    for cc in nx.connected_components(graph):
        model_cc = cc.intersection(model_genes)
        if len(model_cc) > 0:
            if model_cc == model_genes:
                if len(model_cc) == 1 or has_short_path(graph, list(model_cc), unit.functional_unit.transitivity+1):
                    counter[0] += 1
                    sys_state_in_org = "strict"
                else:
                    counter[1] += 1
                    sys_state_in_org = "extended"
            else:
                counter[2] += 1
                sys_state_in_org = "split"
            for cc_gene in cc:
                fam_annot = ""
                fam_sec = []
                if cc_gene.family in gene_fam2mod_fam:
                    metasource, metaid = unit.get_metainfo(cc_gene.family)
                    if metaid != 0:
                        for mod_family in gene_fam2mod_fam[cc_gene.family]:
                            avail_name = {mod_family.name}.union(mod_family.exchangeable)
                            metadata = cc_gene.family.get_metadata(metasource, metaid)
                            if metadata.protein_name in avail_name:
                                fam_annot = metadata.protein_name
                            elif "secondary_name" in metadata.fields:
                                for name in metadata.secondary_name.split(","):
                                    if name in avail_name:
                                        fam_sec.append(name)
                        # print(cc_gene.family.name, cc_gene.family.named_partition)
                        partitions.add(cc_gene.family.named_partition) 
                #     else:
                #         print(f"Warning: No metadata found for family {cc_gene.family.name} in unit {unit.name}")
                # else:
                #     print(f"Warning: No model family found for gene family {cc_gene.family.name} in unit {unit.name}")    
                fam_sec = ",".join(fam_sec)
                projection.append(write_projection_line(cc_gene))
            # print(projection)
            sub_id += 1
    return projection, counter, conciliate_partition(partitions)


def compute_genes_graph(families: Set[GeneFamily], organism: Organism,
                        unit: SystemUnit) -> Tuple[nx.Graph, Set[Gene]]:
    """
    Compute the genes graph for a given genomic context in an organism.

    Args:
        families (Set[GeneFamily]): Set of gene families.
        organism (Organism): The organism of interest.
        unit (SystemUnit): The unit of interest.

    Returns:
        nx.Graph: A genomic context graph for the given organism.
    """
    genes_graph = nx.Graph()
    for family in families:
        genes_graph.add_nodes_from({gene for gene in family.genes if gene.organism == organism})
    mod_fam = set(unit.models_families)
    model_genes = set()
    for gene in sorted(genes_graph.nodes, key=lambda x: x.position):
        if gene.family in mod_fam:
            model_genes.add(gene)
        if gene.position < gene.contig.number_of_genes:
            # print(unit.functional_unit.window, gene.position, gene.contig.number_of_genes)
            right_genes = gene.contig.get_genes(begin=gene.position,
                                                end=gene.position + unit.functional_unit.window + 1,
                                                outrange_ok=True)
        else:
            right_genes = [gene]

        left_genes = gene.contig.get_genes(begin=gene.position - unit.functional_unit.window, end=gene.position + 1,
                                           outrange_ok=True)
        # print(gene.position)
        # print(f"Left genes: {[g.position for g in left_genes]}")
        # print(f"Right genes: {[g.position for g in right_genes]}")
        for l_idx, l_gene in enumerate(left_genes, start=1):
            # print(f"Left index: {l_idx}, window: {unit.functional_unit.window}")
            # if l_gene in genes_graph.nodes:
            for t_gene in left_genes[l_idx:]:
                # print(f"Left gene: {l_gene.position}, Target gene: {t_gene.position}")
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    if t_gene.strand == l_gene.strand:
                        genes_graph.add_edge(t_gene, l_gene, transitivity=l_idx)
                else:
                    genes_graph.add_edge(t_gene, l_gene, transitivity=l_idx)

        for r_idx, r_gene in enumerate(right_genes, start=1):
            # print(f"Right index: {r_idx}, window: {unit.functional_unit.window}")
            # if r_gene in genes_graph.nodes:
            for t_gene in right_genes[r_idx:]:
                # print(f"Right gene: {r_gene.position}, Target gene: {t_gene.position}") 
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    if t_gene.strand == r_gene.strand:
                        genes_graph.add_edge(t_gene, r_gene, transitivity=r_idx)
                else:
                    genes_graph.add_edge(t_gene, r_gene, transitivity=r_idx)
    # print("Before",len(genes_graph.nodes))
    genes_graph.remove_nodes_from([n for n in genes_graph.nodes if n.family not in unit.families])
    # print([n for n in genes_graph.nodes if n.family not in unit.families])
    # print("After",len(genes_graph.nodes))
    # export_graph_for_gephi(genes_graph, f"synced/gene_graphs/genes_graph_{organism.name}_{unit}.gexf")
    # exit(0)
    # import matplotlib.pyplot as plt
    # # colore noes not in model families
    # for gene in genes_graph.nodes:
    #     if gene not in model_genes:
    #         genes_graph.nodes[gene]['color'] = 'blue'
    #     else:
    #         genes_graph.nodes[gene]['color'] = 'red'
    # nx.draw(genes_graph, node_size=10, with_labels=True, font_size=4, node_color=[genes_graph.nodes[n]['color'] for n in genes_graph.nodes])
    # plt.savefig(f"synced/gene_graphs/genes_graph_{organism.name}_{unit.ID}.png", dpi=300)
    # plt.close()
    # remove nodes which families are not in the model families
    return genes_graph, model_genes
# from ppanggolin.utils import extract_contig_window
# from ppanggolin.context.searchGeneContext import add_edges_to_context_graph

# def compute_genes_graph(families: Set[GeneFamily], organism: Organism,
#                         unit: SystemUnit) -> Tuple[nx.Graph, Set[Gene]]:
    
#     genes_graph = nx.Graph()
#     model_genes = set()
#     for family in families:
#         genes_graph.add_nodes_from({gene for gene in family.genes if gene.organism == organism})
#     mod_fam = set(unit.models_families)
#     for gene in sorted(genes_graph.nodes, key=lambda x: x.position):
#         if gene.family in mod_fam:
#             model_genes.add(gene)
#     for contig in organism.contigs:
#         contig_windows = extract_contig_window(contig.number_of_genes, [g.position for g in genes_graph.nodes if g.contig == contig],
#                 window_size=unit.functional_unit.window,is_circular=contig.is_circular)
#         add_edges_to_context_graph(genes_graph, contig, contig_windows, unit.functional_unit.transitivity)
        
#     return genes_graph, model_genes

def export_graph_for_gephi(graph: nx.Graph, filename: str):
    """
    Export a NetworkX graph to GEXF format for Gephi, removing problematic attributes.
    """
    G = graph.copy()
    # Remove non-serializable node attributes
    for n, attrs in G.nodes(data=True):
        for k in list(attrs.keys()):
            v = attrs[k]
            if not isinstance(v, (str, int, float, bool)):
                del G.nodes[n][k]
    # Remove non-serializable edge attributes
    for u, v, attrs in G.edges(data=True):
        for k in list(attrs.keys()):
            val = attrs[k]
            if not isinstance(val, (str, int, float, bool)):
                del G.edges[u, v][k]
    nx.write_gexf(G, filename)
    print(f"Graph exported for Gephi: {filename}")


def find_gene_location(gene, nested_dict):
    for category, gf_dict in nested_dict.items():
        if gene in [gf_id for gf_id, meta in gf_dict.items()]:
            return category
    return None

def unit_projection(unit: SystemUnit, gf2fam: Dict[GeneFamily, set[Family]], fam2source: Dict[str, str],
                    fam_index: Dict[GeneFamily, int], association: List[str] = None
                    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system unit onto all organisms in a pangenome.

    Args:
        unit (SystemUnit): The system unit to project.
        gf2fam (Dict[str, set[Family]]): Dictionary linking a pangenome gene family to a model family.
        fam2source (Dict[str, str]): Dictionary linking a model family to his source.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    pangenome_projection, organisms_projection = [], []
    matrix, m, _ = get_gfs_matrix_combination(set(unit.models_families), gf2fam, fam2source)
    for organism in unit.models_organisms:
        org_fam = {fam for fam in unit.families if organism.bitarray[fam_index[fam]] == 1}
        org_mod_fam = org_fam & set(unit.models_families)
        # print(set(unit.families), set(unit.models_families))
        # print(org_fam - org_mod_fam)
        # exit(0)
        # print([find_gene_location(gene, m) for gene in org_mod_fam])
        filtered_matrix = matrix[list({gf.name for gf in org_mod_fam})]
        # print(filtered_matrix)
        # exit(0) 
        if check_needed_families(filtered_matrix, unit.functional_unit):
            pan_proj = [unit.name, organism.name, ",".join(sorted([x.name for x in org_mod_fam])),
                        ",".join(sorted([x.name for x in org_fam - org_mod_fam]))]
            genes_graph, model_genes = compute_genes_graph(org_mod_fam, organism, unit)
            org_proj, counter, partition = project_unit_on_organisms(genes_graph, unit, organism, model_genes,
                                                                     gf2fam, association)
            pangenome_projection.append(pan_proj + [partition, len(org_fam) / len(unit)] + counter)
            # print(len(unit)==len(set(unit.families)))
            # print(pan_proj + [partition, len(org_fam) / len(unit)] + counter)
            # if 'RGPs' in association:
            #     rgps = {rgp.name for rgp in unit.regions if rgp.organism == organism}
            #     if len(rgps) == 1:
            #         pangenome_projection[-1].extend(rgps)
            #     elif len(rgps) > 1:
            #         join_rgps = [','.join(rgps)]
            #         pangenome_projection[-1].extend(join_rgps)
            # if 'spots' in association:
            #     spots = {str(spot) for spot in unit.spots if organism in spot.organisms}
            #     if len(spots) == 1:
            #         pangenome_projection[-1].extend(spots)
            #     elif len(spots) > 1:
            #         join_spots = [','.join(spots)]
            #         pangenome_projection[-1].extend(join_spots)
            organisms_projection += org_proj
        # else:
        #     print(f"Warning: No needed families found for unit {unit.name} in organism {organism.name}. ")
    return pd.DataFrame(pangenome_projection).drop_duplicates(), pd.DataFrame(organisms_projection).drop_duplicates()


def system_projection(system: System, fam_index: Dict[GeneFamily, int], gene_families: Set[GeneFamily],
                      gene_family2family: Dict[GeneFamily, Set[Family]], fam2source: Dict[str, str],
                      association: List[str] = None
                      ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system onto all organisms in a pangenome.

    Args:
        system (System): The system to project.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        gene_families (Set[GeneFamily]): The set of gene families that code for the model corresponding to system.
        gene_family2family (Dict[GeneFamily, Set[Family]]): Dictionary linking a gene family to model families.
        fam2source (Dict[str, str]): Dictionary linking a model family to his source.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for systems: {system.name}")
    begin = time.time()
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    gfs = gene_families & set(system.families)
    # print("check", gfs==set(system.families), gfs, set(system.families))
    gf2fam = {gf: fam for gf, fam in gene_family2family.items() if gf in gfs}
    for unit in system.units:
        # print(f"Projecting unit: {unit} for system: {system}")
        unit_pan_proj, unit_org_proj = unit_projection(unit, gf2fam, fam2source, fam_index, association)
        pangenome_projection = pd.concat([pangenome_projection, unit_pan_proj], ignore_index=True)
        organisms_projection = pd.concat([organisms_projection, unit_org_proj], ignore_index=True)
    # print(system.ID)
    if len(system) == 1:
        pangenome_projection = pd.concat([pd.DataFrame([[system.ID] * pangenome_projection.shape[0],
                                                        [system.name] * pangenome_projection.shape[0]]).T,
                                          pangenome_projection], axis=1, ignore_index=True)
        organisms_projection = pd.concat([pd.DataFrame([[system.ID] * organisms_projection.shape[0],
                                                        [system.name] * organisms_projection.shape[0]]).T,
                                          organisms_projection], axis=1, ignore_index=True)
        print(f"System {system.ID} Done")

        # if int(system.ID) == 1:
            # write to file for debugging
            # pangenome_projection.to_csv("system_1_pan.tsv", sep="\t", index=False)
            # organisms_projection.to_csv("system_1_org.tsv", sep="\t", index=False)
    else:
        print(f"WARNING: System {system.ID} has multiple units, skipping ID and name projection")

    logging.getLogger("PANORAMA").debug(f"System projection done for {system.name} in {time.time() - begin} seconds")
    return pangenome_projection.drop_duplicates(), organisms_projection.drop_duplicates()


def project_pangenome_systems(pangenome: Pangenome, system_source: str, fam_index: Dict[GeneFamily, int],
                              association: List[str] = None, canonical: bool = False, threads: int = 1,
                              lock: Lock = None, disable_bar: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project systems onto all organisms in a pangenome.

    Args:
        pangenome (Pangenome): The pangenome to project.
        system_source (str): Source of the systems to project.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').
        threads (int, optional): Number of threads available (default is 1).
        lock (Lock, optional): Global lock for multiprocessing execution (default is None).
        disable_bar (bool, optional): Disable progress bar (default is False).

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projections for each organism and the pangenome.
    """
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    meta2fam = get_metadata_to_families(pangenome, pangenome.systems_sources_to_metadata_source()[system_source])
    sys2fam_context = {}

    for system in pangenome.systems:  # Search association now to don't repeat for same model and different system
        if system.model.name not in sys2fam_context:
            gene_families, gf2fam, fam2source = dict_families_context(system.model, meta2fam)
            sys2fam_context[system.model.name] = (gene_families, gf2fam, fam2source)
        if canonical:
            for canonic in system.canonical:
                if canonic.model.name not in sys2fam_context:
                    gene_families, gf2fam, fam2source = dict_families_context(canonic.model, meta2fam)
                    sys2fam_context[canonic.model.name] = (gene_families, gf2fam, fam2source)

    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        logging.getLogger("PANORAMA").info(f'Begin system projection for source : {system_source}')
        with tqdm(total=pangenome.number_of_systems(system_source, with_canonical=canonical), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            print([sys.ID for sys in pangenome.get_system_by_source(system_source)])
            for sys in pangenome.get_system_by_source(system_source):
                for unit in sys.units:
                    print({u for u in unit.families} - {u for u in unit.models_families})
            # print([f for unit in pangenome.get_system_by_source(system_source) for f in unit.models_families])
            # write file with column1 unit and column2 all families in the unit
            # with open("hash_units_and_families.tsv", "w") as out_file:
            #     out_file.write("Unit\tFamilies\n")  # Header (optional)

            #     for unit in pangenome.get_system_by_source(system_source):
            #         families = ",".join(map(str, [f.name for f in unit.models_families]))
            #         out_file.write(f"{unit.name}\t{families}\n")

            # exit(0)
            for system in pangenome.get_system_by_source(system_source):
                # print(f"Projecting system: {system.name}", system)
                gene_families, gf2fam, fam2source = sys2fam_context[system.model.name]
                future = executor.submit(system_projection, system, fam_index, gene_families,
                                         gf2fam, fam2source, association)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                if canonical:
                    for canonic in system.canonical:
                        gene_families, gf2fam, fam2source = sys2fam_context[canonic.model.name]
                        future = executor.submit(system_projection, system, fam_index, gene_families,
                                                 gf2fam, fam2source, association)
                        future.add_done_callback(lambda p: progress.update())
                        futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat([pangenome_projection, result[0]], ignore_index=True)
                organisms_projection = pd.concat([organisms_projection, result[1]], ignore_index=True)
    pan_cols_name = ["system number", "system name", "functional unit name", "organism", "model_GF", "context_GF",
                     "partition", "completeness", "strict", "extended", "split"]
    org_cols_name = ["system number", "system name", "functional unit name", "subsystem number", "organism",
                     "gene family", "partition", "annotation", "secondary_names", "gene.ID", "gene.name", "contig",
                     "start", "stop", "strand", "is_fragment", "genomic organization", "product"]
    # if 'RGPs' in association:
    #     pan_cols_name += ['RGPs']
    #     org_cols_name += ['RGPs']
    # if 'spots' in association:
    #     pan_cols_name += ['spots']
    #     org_cols_name += ['spots']
    pangenome_projection.columns = pan_cols_name
    pangenome_projection.sort_values(
        by=["system number", "system name", "functional unit name", "organism", "completeness"],
        ascending=[True, True, True, True, True],
        inplace=True)  # TODO Try to order system number numerically
    organisms_projection.columns = org_cols_name
    organisms_projection.sort_values(
        by=["system name", "system number", "subsystem number", "functional unit name",
            "organism", "contig", "start", "stop"],
        ascending=[True, True, False, True, True, True, True, True],
        inplace=True)
    logging.getLogger("PANORAMA").debug('System projection done')
    # print(pangenome_projection.head())
    return pangenome_projection, organisms_projection


def _custom_agg(series: pd.Series, unique: bool = False):
    """
    Aggregate a column

    Args:
        series: series to aggregate
        unique: whether to return unique values or not

    Returns:
        The aggregated series
    """
    if unique:
        values = [set(x) for x in series.replace('', pd.NA).dropna().str.split(',')]
        # print(values)
        if len(values) == 0:
            return ''
        else:  # len(values) >1
            return ', '.join(sorted(set().union(*values)))
    else:
        values = list(itertools.chain(*[x for x in series.replace('', pd.NA).dropna().str.split(',')]))
        if len(values) == 0:
            return ''
        else:  # len(values) >=1
            return ', '.join(sorted(values))


def custom_agg(series: pd.Series):
    """
    Aggregate a column

    Args:
        series: series to aggregate

    Returns:
        The aggregated series
    """
    return _custom_agg(series, unique=False)


def custom_agg_unique(series: pd.Series):
    """
    Aggregate a column

    Args:
        series: series to aggregate

    Returns:
        The aggregated series
    """
    # print(series)
    return _custom_agg(series, unique=True)


def get_partition(series: pd.Series):
    """

    Args:
        series:

    Returns:

    """

    final_partition = ""
    for partition in series.unique().tolist():
        if partition != final_partition:
            if final_partition == "":
                final_partition = partition
            else:
                if final_partition == "persistent":
                    if partition == "shell":
                        final_partition = "persistent/shell"
                    elif partition == "cloud":
                        final_partition = "persistent/cloud"
                    elif partition == "accessory":
                        final_partition = "persistent/accessory"
                elif final_partition == "shell":
                    if partition == "persistent":
                        final_partition = "persistent/shell"
                    elif partition == "cloud":
                        final_partition = "accessory"
                    elif partition == "persistent/cloud":
                        final_partition = "persistent/accessory"
                elif final_partition == "cloud":
                    if partition == "persistent":
                        final_partition = "persistent/cloud"
                    elif partition == "shell":
                        final_partition = "accessory"
                    elif partition == "persistent/shell":
                        final_partition = "persistent/accessory"
                elif final_partition == "accessory":
                    if partition == "persistent":
                        final_partition = "persistent/accessory"
                elif final_partition == "persistent/shell":
                    if partition == "cloud":
                        final_partition = "persistent/accessory"
                elif final_partition == "persistent/cloud":
                    if partition == "shell":
                        final_partition = "persistent/accessory"
    return final_partition


def extract_numeric_for_sorting(val) -> float:
    """
    Function to extract the numeric value for sorting while keeping the original value

    Args:
        val: the value

    Returns:
        float: the numeric value
    """
    try:
        # If it's a simple number, return it for sorting
        return float(val)
    except ValueError:
        # If it's a list of numbers separated by commas, return the smallest number for sorting
        if ',' in val:
            parts = [float(x) for x in val.replace('"', '').split(',')]
            return min(parts)  # Take the minimum for sorting
        return float('inf')  # If it cannot be converted, place it at the end


def get_org_df(org_df: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    """
    Get the reformated projection dataframe for an organism

    Args:
        org_df: Dataframe for the corresponding organism

    Returns:
        pd.DataFrame: Dataframe reformated for an organism
    """
    org_name = org_df["organism"].unique()[0]
    org_df = org_df.drop(columns=["organism"])
    org_df_cols = org_df.columns.tolist()

    # print(org_df)
    # Create a temporary column for sorting based on the numeric values extracted
    org_df["sort_key"] = org_df["system number"].apply(extract_numeric_for_sorting)

    # Sort the DataFrame using the temporary column, but keep the original values
    org_df_sorted = org_df.sort_values(by=["sort_key", "system name", "start", "stop"],
                                       ascending=[True, True, True, True]).drop(columns=["sort_key"])
    org_df_grouped = org_df_sorted.groupby(["gene family", "system name", "functional unit name", "gene.ID", "start"],
                                           as_index=False)
    
    agg_dict = {"system number": custom_agg_unique, "subsystem number": custom_agg, "partition": custom_agg_unique,
                "annotation": custom_agg, "secondary_names": custom_agg, "contig": custom_agg_unique,
                "gene.name": custom_agg_unique, "stop": custom_agg_unique, "strand": custom_agg_unique,
                "is_fragment": custom_agg_unique, "genomic organization": custom_agg, "product": custom_agg}
    if "RGPs" in org_df_cols:
        agg_dict["RGPs"] = custom_agg_unique
    if "spots" in org_df_cols:
        agg_dict["spots"] = custom_agg_unique
    org_df_grouped = org_df_grouped.agg(agg_dict)
    org_df_grouped = org_df_grouped[org_df_cols]
    # print(org_df_grouped)
    # Create a temporary column for sorting based on the numeric values extracted
    org_df_grouped["sort_key"] = org_df_grouped["system number"].apply(extract_numeric_for_sorting)

    # Sort the DataFrame using the temporary column, but keep the original values
    org_df_grouped_sorted = org_df_grouped.sort_values(by=["sort_key", "system name", "start", "stop"],
                                                       ascending=[True, True, True, True]).drop(columns=["sort_key"])

    return org_df_grouped_sorted, org_name


def write_projection_systems(output: Path, pangenome_projection: pd.DataFrame, organisms_projection: pd.DataFrame,
                             organisms: List[str] = None, threads: int = 1, force: bool = False, disable_bar: bool = False):
    """
    Write the projected systems to output files.

    Args:
        output (Path): Path to the output directory.
        pangenome_projection (pd.DataFrame): DataFrame containing the pangenome projection.
        organisms_projection (pd.DataFrame): DataFrame containing the organism projections.
        organisms (List[str], optional): List of organisms to project (default is all organisms).
        force (bool, optional): Force write to the output directory (default is False).

    Returns:
        None
    """

    proj_dir = mkdir(output / "projection", force=force)
    if organisms is not None:
        pangenome_projection = pangenome_projection[~pangenome_projection["organism"].isin(organisms)]
        organisms_projection = organisms_projection[~organisms_projection["organism"].isin(organisms)]

    with ProcessPoolExecutor(max_workers=threads, mp_context=get_context("fork")) as executor:
        futures = []
        for organism_name in pangenome_projection["organism"].unique():
            org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
            future = executor.submit(get_org_df, org_df)
            futures.append(future)

        for future in tqdm(as_completed(futures), total=len(pangenome_projection["organism"].unique()),
                           unit='organisms', disable=disable_bar, desc="System projection on organisms"):
            try:
                organism_df, organim_name = future.result()
                print(f"Writing projection for organism: {organim_name}")
                organism_df.to_csv(proj_dir / f"{organim_name}.tsv", sep="\t", index=False)
            except Exception as e:
                print(f"[ERROR] A future failed: {e}")

    pan_df_col = pangenome_projection.columns.tolist()
    pangenome_grouped = pangenome_projection.groupby(by=["system number", "system name"], as_index=False)
    agg_dict = {"functional unit name": custom_agg_unique, "organism": custom_agg_unique,
                "model_GF": custom_agg_unique, "context_GF": custom_agg_unique,
                "partition": get_partition, "completeness": 'mean',
                "strict": 'sum', "extended": 'sum', "split": 'sum'}
    if "RGPs" in pan_df_col:
        agg_dict["RGPs"] = custom_agg_unique
    if "spots" in pan_df_col:
        agg_dict["spots"] = custom_agg_unique

    pangenome_grouped = pangenome_grouped.agg(agg_dict)
    pangenome_grouped = pangenome_grouped[pan_df_col]

    # Create a temporary column for sorting based on the numeric values extracted
    pangenome_grouped["sort_key"] = pangenome_grouped["system number"].apply(extract_numeric_for_sorting)

    # Sort the DataFrame using the temporary column, but keep the original values
    pangenome_sorted = pangenome_grouped.sort_values(by=["sort_key", "system name", "organism"],
                                                     ascending=[True, True, True]).drop(columns=["sort_key"])

    pangenome_sorted.to_csv(output / 'systems.tsv', sep="\t", index=False)
