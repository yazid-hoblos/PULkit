#!/usr/bin/env python3

"""
This module provides utility functions to detect and write biological systems in pangenomes.
"""

# default libraries
from __future__ import annotations
import logging
from typing import Dict, Iterable, List, Set, Tuple
from collections import defaultdict

# installed libraries
from itertools import combinations
import numpy as np
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism
from ppanggolin.metadata import Metadata

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.pangenomes import Pangenome

pd.options.mode.copy_on_write = True  # Remove when pandas3.0 available. See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy


def filter_global_context(graph: nx.Graph, jaccard_threshold: float = 0.8) -> nx.Graph[GeneFamily]:

    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = len(data["genomes"]) / len(set(f1.organisms))
        f2_proportion = len(data["genomes"]) / len(set(f2.organisms))

        # Update the local copy of the edge data
        data.update({
            'f1': f1.name,
            'f2': f2.name,
            'f1_jaccard_gene': f1_proportion,
            'f2_jaccard_gene': f2_proportion
        })

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def filter_local_context(graph: nx.Graph, organisms: Set[Organism],
                         jaccard_threshold: float = 0.8) -> nx.Graph[GeneFamily]:
    print("Filtering local context with Jaccard threshold:", jaccard_threshold)
    print("******************************")
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """
    def get_gene_proportion(gene_family: GeneFamily) -> float:
        """Returns the Jaccard gene proportion for a given gene family."""
        # Compute the proportion if not cached
        gf_orgs = set(gene_family.organisms).intersection(organisms)
        if len(gf_orgs) > 0:
            return len(data["genomes"].intersection(organisms)) / len(gf_orgs)
        else:
            return 0

    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = get_gene_proportion(f1)
        f2_proportion = get_gene_proportion(f2)

        # Update the local copy of the edge data
        data.update({
            'f1': f1.name,
            'f2': f2.name,
            'f1_jaccard_gene': f1_proportion,
            'f2_jaccard_gene': f2_proportion
        })

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def get_number_of_mod_fam(gene_family: GeneFamily, gene_fam2mod_fam: Dict[GeneFamily, Set[Family]]) -> int:
    """Gets the number of model families associated with the gene family.

    Args:
        gene_family (GeneFamily): The gene family of interest.

    Returns:
        int: The number of model families associated with the gene family.
    """
    model_families = gene_fam2mod_fam.get(gene_family.name)
    if model_families is not None:
        return len(model_families)
    else:
        return 0


def check_for_families(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                       mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit
                       ) -> Tuple[bool, Dict[GeneFamily, Tuple[str, int]]]:
    """
    Checks if there are forbidden conditions in the families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.
        func_unit (FuncUnit): Functional unit to check against.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    mandatory_list = list(map(lambda x: x.name, func_unit.mandatory))
    accessory_list = list(map(lambda x: x.name, func_unit.accessory))
    forbidden_list = list(map(lambda x: x.name, func_unit.forbidden))

    gf2fam2meta_info = defaultdict(dict)

    for gf in sorted(gene_families, key=lambda n: get_number_of_mod_fam(n, gene_fam2mod_fam)):
        for family in gene_fam2mod_fam[gf]:
            avail_name = {family.name}.union(family.exchangeable)
            if ((family.presence == 'mandatory' and family.name in mandatory_list) or
                    (family.presence == 'accessory' and family.name in accessory_list) or
                    (family.presence == 'forbidden' and family.name in forbidden_list)):
                for meta_id, metadata in gf.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
                    if metadata.protein_name in avail_name:
                        gf2fam2meta_info[gf][family] = (mod_fam2meta_source[family.name], meta_id, metadata.score)
                    elif "secondary_name" in metadata.fields:
                        if any(name in avail_name for name in metadata.secondary_name.split(",")):
                            gf2fam2meta_info[gf][family] = (mod_fam2meta_source[family.name], meta_id, metadata.score)

    gf2meta_info = {}
    mandatory_seen = set()
    accessory_seen = set()
    found_forbidden = False
    for gf, fam2meta_info in sorted(gf2fam2meta_info.items(), key=lambda x: len(x[1])):
        sorted_fam2meta_info = list(sorted(fam2meta_info.items(), key=lambda x: x[1][2], reverse=True))
        add = False
        for family, meta_info in sorted_fam2meta_info:
            if family.presence == "mandatory" and family.name not in mandatory_seen:
                mandatory_seen.add(family.name)
                gf2meta_info[gf] = meta_info[:-1]
                add = True
                break
            elif family.presence == "accessory" and family.name not in accessory_seen:
                accessory_seen.add(family.name)
                gf2meta_info[gf] = meta_info[:-1]
                add = True
                break
            elif family.presence == "forbidden" and family.name in forbidden_list:
                found_forbidden = True
                break  # No need to continue if a forbidden family is found

        if found_forbidden:
            break

        if not add:
            # If nothing was added and no forbidden family was found, default to first family
            _, meta_info = sorted_fam2meta_info[0]
            gf2meta_info[gf] = meta_info[:-1]

    if found_forbidden:
        return False, {}
    else:
        if ((len(mandatory_seen) >= func_unit.min_mandatory or func_unit.min_mandatory == -1) and
                (len(mandatory_seen | accessory_seen) >= func_unit.min_total or func_unit.min_total == -1)):
            return True, gf2meta_info
        else:
            return False, {}


def get_gfs_matrix_combination(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                               mod_fam2meta_source: Dict[str, str]
                               ) -> Tuple[pd.DataFrame, Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                                          Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]]:
    """
    Build a matrix of association between gene families and families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]): Dictionary linking gene families to model families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking families to metadata sources.

    Returns:
        pd.DataFrame: Matrix of association between gene families and families.
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to mandatory families
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to accessory families
    """

    def add_metadata_to_dict(presence_gfs2metadata):
        """
        Associate gene families to metadata in a dictionary

        Args:
            presence_gfs2metadata: Dictionary to save gene families and their metadata
        """
        # print(gene_family.get_metadata_by_source(mod_fam2meta_source[family.name]).items())
        for meta_id, metadata in gene_family.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
            if metadata.protein_name in avail_name:
                if gene_family in presence_gfs2metadata[family.name]:
                    _, current_metadata = presence_gfs2metadata[family.name][gene_family]
                    if metadata.score > current_metadata.score:
                        presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
                else:
                    presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
            elif "secondary_name" in metadata.fields:
                if any(name in avail_name for name in metadata.secondary_name.split(",")):
                    if gene_family in presence_gfs2metadata[family.name]:
                        _, current_metadata = presence_gfs2metadata[family.name][gene_family]
                        if metadata.score > current_metadata.score:
                            presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
                    else:
                        presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)

    mandatory_gfs2metadata = defaultdict(dict)
    accessory_gfs2metadata = defaultdict(dict)
    gfs = set()
    for gene_family in gene_families:
        for family in gene_fam2mod_fam[gene_family]:
            avail_name = {family.name}.union(family.exchangeable)
            if family.presence == "mandatory":
                add_metadata_to_dict(mandatory_gfs2metadata)
                gfs.add(gene_family)
            elif family.presence == "accessory":
                add_metadata_to_dict(accessory_gfs2metadata)
                gfs.add(gene_family)
    # print(mandatory_gfs2metadata)
    # print(accessory_gfs2metadata)
    fams = list(mandatory_gfs2metadata.keys()) + list(accessory_gfs2metadata.keys())
    score_matrix = np.zeros((len(fams), len(gfs)))
    gfs = list(gfs)
    for i, fam in enumerate(fams):
        if fam in mandatory_gfs2metadata:
            gfs2metadata = mandatory_gfs2metadata[fam]
        else:
            gfs2metadata = accessory_gfs2metadata[fam]
        for j, gf in enumerate(gfs):
            score_matrix[i, j] = 1 if gf in gfs2metadata else 0
    return (pd.DataFrame(score_matrix, index=fams, columns=[gf.name for gf in gfs]),
            mandatory_gfs2metadata, accessory_gfs2metadata)


def bitset_from_row(row) -> int:
    """
    Converts a row of the binary matrix into an integer bitset.
    Each bit represents the presence of an individual for that characteristic.
    Args:
        row: A row of the binary matrix

    Returns:
        int: integer bitset representing the presence of an individual for that characteristic
    """
    bitset = 0
    for i, val in enumerate(row):
        if val == 1:
            bitset |= (1 << i)  # Met le i-ème bit à 1
    return bitset


def search_comb(comb: Tuple[int, ...], mandatory_bitsets: List[int],
                accessory_bitsets: List[int]) -> Tuple[Dict[int, int], int, int]:
    """
    Search for a working combination in bitsets

    Args:
        comb (Tuple[int, ...]): Combination of gene families index
        mandatory_bitsets (List[int]): Bitsets with index of mandatory families
        accessory_bitsets (List[int]): Bitsets with index of accessory families

    Returns:
        Dict[int, int]: Association between gene families and families
        Integer: Number of mandatory families covered
        Integer: Number of accessory families covered
    """
    covered_families = set()
    gf2fam = {}

    covered_mandatory = 0
    for family_index, bitset in enumerate(mandatory_bitsets):
        covering_gfs = set(ind for ind in comb if (1 << ind) & bitset)
        if len(covering_gfs) >= 1:
            for gf in covering_gfs:
                gf2fam[gf] = family_index
            if family_index not in covered_families:
                covered_mandatory += 1
                covered_families.add(family_index)

    covered_accessory = 0
    for family_index, bitset in enumerate(accessory_bitsets, start=len(mandatory_bitsets)):
        covering_gfs = set(ind for ind in comb if (1 << ind) & bitset)
        if len(covering_gfs) >= 1:
            for gf in covering_gfs:
                gf2fam[gf] = family_index
            if family_index not in covered_families:
                covered_accessory += 1
                covered_families.add(family_index)

    return gf2fam, covered_mandatory, covered_accessory


def find_combinations(matrix: pd.DataFrame, func_unit: FuncUnit) -> List[Tuple[Set[str], Dict[str, str]]]:
    """
    Search working combination of gene families that respect families presence absence model rules

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        List: List of working combination
            Set[str]: Set of selected gene families name that correspond to a working combination
            Dict[str, str]: Association between gene families and families for the combination.
    """

    mandatory = {fam.name for fam in func_unit.mandatory}
    mandatory_indices = [i for i, fam in enumerate(matrix.index.values) if fam in mandatory]

    if len(mandatory_indices) < func_unit.min_mandatory or matrix.shape[0] < func_unit.min_total:
        return []

    bitsets = [bitset_from_row(matrix.iloc[i, :]) for i in range(matrix.shape[0])]
    mandatory_bitsets = [bitsets[i] for i in mandatory_indices]
    accessory_bitsets = [bitsets[i] for i in range(matrix.shape[0]) if i not in mandatory_indices]

    solutions = []
    for size in sorted(range(func_unit.min_total, matrix.shape[1] + 1), reverse=True):
        for comb in combinations(range(matrix.shape[1]), size):
            gf2fam, covered_mandatory, covered_accessory = search_comb(comb, mandatory_bitsets, accessory_bitsets)
            if (covered_mandatory >= func_unit.min_mandatory and
                    covered_mandatory + covered_accessory >= func_unit.min_total):
                solutions.append(({matrix.columns[i] for i in comb},
                                  {matrix.columns[j]: matrix.index[i] for j, i in gf2fam.items()}))

    return solutions


def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Search if it exists a combination of gene families that respect families presence absence unit rules

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if it exists, False otherwise
    """
    # logging.getLogger("PANORAMA").debug("Checking if the functional unit can be satisfied")
    def is_subset(comb1: Tuple[int, ...], comb2: Tuple[int, ...]) -> bool:
        """
        Check if a combination of gene families is subset of another one

        Args:
            comb1 (Tuple[int, ...]): Combination to check if it's a subset
            comb2 (Tuple[int, ...]): Combination supposed bigger

        Returns:
            Boolean: True if first combination is a subset of the second one, False otherwise
        """
        return set(comb1).issubset(comb2)

    mandatory = {fam.name for fam in func_unit.mandatory}
    mandatory_indices = [i for i, fam in enumerate(matrix.index.values) if fam in mandatory]

    # if (matrix.iloc[0] == 0).all() or (matrix.iloc[1] == 0).all():
    #     print("Matrix is empty")

    if len(mandatory_indices) < func_unit.min_mandatory or matrix.shape[0] < func_unit.min_total:
        return False

    # if (matrix == 0).all(axis=1).any():
        # return False

    bitsets = [bitset_from_row(matrix.iloc[i, :]) for i in range(matrix.shape[0])]
    mandatory_bitsets = [bitsets[i] for i in mandatory_indices]
    accessory_bitsets = [bitsets[i] for i in range(matrix.shape[0]) if i not in mandatory_indices]

    _, covered_mandatory, covered_accessory = search_comb(tuple(range(matrix.shape[1])),mandatory_bitsets,accessory_bitsets)

    # print(covered_mandatory >= func_unit.min_mandatory and
    #         covered_mandatory + covered_accessory >= func_unit.min_total)
    return (covered_mandatory >= func_unit.min_mandatory and
            covered_mandatory + covered_accessory >= func_unit.min_total)
    #     return True
    # return False

    not_valid = []
    for size in sorted(range(func_unit.min_total, matrix.shape[1] + 1), reverse=True):
        num_combs = sum(1 for _ in combinations(range(matrix.shape[1]), size))
        logging.getLogger("PANORAMA").debug(f"Checking {num_combs} combinations of size {size}")
        for comb in combinations(range(matrix.shape[1]), size):
            # If a larger comb can't respect presence absence rules, so a tinier can not too.
            if not any(is_subset(comb, larger_comb) for larger_comb in not_valid):
                _, covered_mandatory, covered_accessory = search_comb(comb, mandatory_bitsets,
                                                                      accessory_bitsets)
                print(covered_mandatory, covered_accessory)
                print(func_unit.min_mandatory, func_unit.min_total)

                if (covered_mandatory >= func_unit.min_mandatory and
                        covered_mandatory + covered_accessory >= func_unit.min_total):
                    print("Found a working combination")
                    return True
                else:
                    not_valid.append(comb)
            # else:
                # print("Still checking for a working combination")
    return False


def get_metadata_to_families(pangenome: Pangenome, sources: Iterable[str]) -> Dict[str, Dict[str, Set[GeneFamily]]]:
    """
    Retrieves a mapping of metadata to sets of gene families for each metadata source.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        sources (iterable of str): List of metadata source names.

    Returns:
        dict: A dictionary where each metadata source maps to another dictionary of metadata to sets of gene families.
    """
    meta2fam = {source: defaultdict(set) for source in sources}
    for source in sources:
        for gf in pangenome.gene_families:
            metadata = gf.get_metadata_by_source(source)
            if metadata is not None:
                for meta in metadata.values():
                    meta2fam[source][meta.protein_name].add(gf)
                    if "secondary_name" in meta.fields and meta.secondary_name != "":
                        for secondary_name in meta.secondary_name.split(','):
                            meta2fam[source][secondary_name].add(gf)
    return meta2fam


def dict_families_context(model: Model, annot2fam: Dict[str, Dict[str, Set[GeneFamily]]]) \
        -> Tuple[Set[GeneFamily], Dict[GeneFamily, Set[Family]], Dict[str, str]]:
    """
    Retrieves all gene families associated with the families in the model.

    Args:
        model (Model): Model containing the families.
        annot2fam (dict): Dictionary of annotated families.

    Returns:
        tuple: A tuple containing:
            - Set[GeneFamily]: Gene families of interest in the functional unit.
            - dict: Dictionary linking gene families to their families.
            - dict: Dictionary linking families to their sources.
    """
    gene_families = set()
    gf2fam = defaultdict(set)
    fam2source = {}
    for fam_model in model.families:
        for source, annotation2families in annot2fam.items():
            if fam_model.name in annotation2families:
                for gf in annotation2families[fam_model.name]:
                    gene_families.add(gf)
                    gf2fam[gf].add(fam_model)
                    if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                        logging.getLogger("PANORAMA").warning("Two families have the same protein name for different "
                                                              "sources. First source encountered will be used.")
                    else:
                        fam2source[fam_model.name] = source

        for exchangeable in fam_model.exchangeable:
            for source, annotation2families in annot2fam.items():
                if exchangeable in annotation2families:
                    for gf in annotation2families[exchangeable]:
                        gene_families.add(gf)
                        gf2fam[gf].add(fam_model)
                        if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                            logging.getLogger("PANORAMA").warning(
                                "Two families have the same protein name for different "
                                "sources. First source encountered will be used.")
                        else:
                            fam2source[fam_model.name] = source
    # print families in gf2fam with more than one corresponding gene family
    # for gf, families in gf2fam.items():
    #     if len(families) > 1:
    #         print(f"Gene family {gf.name} is associated with multiple families: {', '.join(f.name for f in families)}")
    return gene_families, gf2fam, fam2source
