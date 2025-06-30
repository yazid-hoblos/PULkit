#!/usr/bin/env python3
# coding:utf-8

"""
This module provides function to translate models from different databases
"""

# default libraries
import re
from random import choice
from string import digits
from typing import Callable, Dict, List, Set, Union
import logging
from pathlib import Path
from lxml import etree as et
import lxml.etree
import yaml
import json

from networkx.algorithms.cluster import transitivity
# installed libraries
from tqdm import tqdm
import pandas as pd
from numpy import nan
from pyhmmer.plan7 import HMM

# local libraries
from panorama.utils import mkdir
from panorama.utility.genInput import create_hmm_list_file, read_hmm, write_hmm, gen_acc


known_sources = ["padloc", "defense-finder", "CONJScan", "TXSScan", "TFFscan"]


def read_yaml(model: Path) -> dict:
    """
    Reads a YAML file and returns its contents as a Python object.

    Args:
        model (str): The path to the YAML file to be read.

    Raises:
        IOError: If there is a problem opening the file.
        Exception: If there is an unexpected problem reading the file.

    Returns:
        dict: The contents of the YAML file, as a Python object.
    """
    try:
        model_file = open(model, "r")
    except IOError as ioerror:
        raise IOError(f"Problem to open {model}. The following error is the direct cause: {ioerror}")
    except Exception as error:
        raise Exception(f"Unexpected problem reading {model}. The following error is the direct cause: {error}")
    else:
        data = yaml.safe_load(model_file)
        model_file.close()
        return data


def read_xml(model: Path) -> et.Element:
    """
    Reads an XML file and returns its contents as a Python object.

    Args:
        model (str): The path to the YAML file to be read.

    Raises:
        IOError: If there is a problem opening the file.
        Exception: If there is an unexpected problem reading the file.

    Returns:
        et.Element: The contents of the YAML file, as a Python object.
    """
    try:
        parser = et.XMLParser(remove_comments=True, resolve_entities=False, no_network=True)
        my_tree = et.parse(model.absolute().as_posix(), parser=parser)
    except IOError as ioerror:
        raise IOError(f"Problem to open {model}. The following error is direct cause\n{ioerror}")
    except Exception as error:
        raise Exception(f"Unexpected Problem to read {model}. The following error is direct cause\n{error}")
    else:
        my_root = my_tree.getroot()
        return my_root


def write_model(path_output: Path, dict_data: dict) -> None:
    """
    This function writes a model to a JSON file.

    Args:
        path_output (Path): The path to the output directory.
        dict_data (dict): The model data to be written to the file.

    Returns:
        None
    """
    with open(f"{path_output}/models/{dict_data['name']}.json", "w") as my_file:
        json.dump(dict_data, my_file, indent=2)


def parse_meta_padloc(meta: Path) -> pd.DataFrame:
    """
    This function parses a PADLOC metadata file and returns a pandas DataFrame containing the metadata.

    Args:
        meta (Path): The path to the PADLOC metadata file.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the metadata.
    """
    def merge_columns(row):
        if pd.isna(row['temp']) and pd.isna(row['secondary.name']):
            return ''
        elif pd.isna(row['temp']):
            return row['secondary.name']
        elif pd.isna(row['secondary.name']):
            return row['temp']
        else:
            return f"{row['temp']},{row['secondary.name']}"

    meta_col_names = ["accession", "name", "protein_name", "secondary_name", "score_threshold",  "eval_threshold",
                      "ieval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
    df = pd.read_csv(meta, sep="\t", usecols=[0, 1, 2, 3, 4, 6, 7, 8], header=0)
    df[['protein_name', 'temp']] = df['protein.name'].str.split('|', expand=True)
    df['secondary_name'] = df.apply(merge_columns, axis=1)
    df = df.drop('protein.name', axis=1)
    df = df.drop('secondary.name', axis=1)
    df = df.drop('temp', axis=1)
    df = df.iloc[:, [0, 1, 6, 7, 3, 4, 5, 2]]
    df.insert(4, 'score_threshold', None)
    df.insert(5, 'eval_threshold', None)
    df.columns = meta_col_names
    df = df.set_index('accession')
    df['description'] = df["description"].fillna('unknown')
    return df


def translate_model_padloc(data_yaml: dict, model_name: str, meta: pd.DataFrame = None, canonical: list = None) -> dict:
    """
    This function translates a PADLOC model from a YAML file into a JSON format that can be used by the PADLOC software.

    Args:
        data_yaml (dict): The PADLOC model data in YAML format.
        model_name (str): The name of the PADLOC model.
        meta (pd.DataFrame, optional): The metadata associated with the model. Defaults to None.
        canonical (list, optional): The canonical genes associated with the model. Defaults to None.

    Returns:
        dict: The translated PADLOC model in JSON format.

    Raises:
        AssertionError: if the meta argument is not a dataframe or if the canonical argument is not a list
        KeyError: If the PADLOC model data contains unexpected keys.
    """
    assert canonical is not None and isinstance(canonical, list)
    assert meta is not None and isinstance(meta, pd.DataFrame)

    def add_families(families_list: List[str], sec_names: List[str], fam_type: str) -> List[Dict[str, str]]:
        """
        Add a list of families to the functional unit

        Args:
            families_list: list of families to add
            sec_names: list of known secondary names for the families
            fam_type: type of families

        Returns:
            List[Dict[str, str]]: list of families filled with information to be a PANORAMA model

        Todo:
            * Try to find an elegant way to manage the cas_adaptation exception
        """
        fam_list = []
        for fam_name in families_list:
            if fam_name not in seen_fam:
                if fam_name == "cas_adaptation":
                    filter_df = meta.loc[meta['secondary_name'] == fam_name]
                    for filter_fam in filter_df['protein_name'].dropna().unique().tolist():
                        if filter_fam not in seen_fam:
                            fam_dict = {'name': filter_fam,
                                        'presence': fam_type}
                            seen_fam.add(filter_fam)
                            if filter_fam in sec_names:
                                filter_df = meta.loc[meta['secondary_name'] == fam_name]
                                fam_dict.update({"exchangeable": filter_df['protein_name'].dropna().unique().tolist()})
                            fam_list.append(fam_dict)
                elif fam_name == "cas_accessory":
                    pass
                elif fam_name != 'NA':
                    fam_dict = {'name': fam_name,
                                'presence': fam_type}
                    seen_fam.add(fam_name)
                    if fam_name in sec_names:
                        filter_df = meta.loc[meta['secondary_name'] == fam_name]
                        fam_dict.update({"exchangeable": filter_df['protein_name'].dropna().unique().tolist()})
                    fam_list.append(fam_dict)
        return fam_list

    padloc_keys = ["maximum_separation", "minimum_core", "minimum_total", "force_strand",
                   "core_genes", "secondary_genes", "neutral_genes", "prohibited_genes"]
    if not all(key in padloc_keys for key in data_yaml.keys()):
        raise KeyError(f"Unexpected key in PADLOC model : {model_name}."
                       f"authorized keys are : {', '.join(padloc_keys)}")
    if "core_genes" not in data_yaml:
        raise KeyError("Core gene must be found in padloc keys")

    data_json = {"name": model_name, 'func_units': [],
                 'parameters': {"transitivity": 0, "window": 1, "min_mandatory": 1, "min_total": 1}
                 }
    if len(canonical) > 0:
        data_json['canonical'] = canonical
    secondary_names = meta['secondary_name'].dropna().unique().tolist()
    func_unit = {'name': model_name, 'presence': 'mandatory', 'same_strand': data_yaml["force_strand"],
                 'parameters':
                     {
                         "transitivity": data_yaml["maximum_separation"] if "maximum_separation" in data_yaml else 0,
                         "min_mandatory": data_yaml["minimum_core"] if "minimum_core" in data_yaml else 1,
                         "min_total": data_yaml["minimum_total"] if "minimum_total" in data_yaml else 1
                     }
                 }
    family_list = list()
    seen_fam = set()
    family_list += add_families(families_list=data_yaml["core_genes"], sec_names=secondary_names,
                                fam_type='mandatory')
    family_list += add_families(families_list=data_yaml["secondary_genes"] if "secondary_genes" in data_yaml else [],
                                sec_names=secondary_names, fam_type='accessory')
    family_list += add_families(families_list=data_yaml["prohibited_genes"] if "prohibited_genes" in data_yaml else [],
                                sec_names=secondary_names, fam_type='forbidden')
    family_list += add_families(families_list=data_yaml["neutral_genes"] if "neutral_genes" in data_yaml else [],
                                sec_names=secondary_names, fam_type='neutral')
    func_unit["families"] = family_list
    data_json['func_units'].append(func_unit)
    return data_json


def search_canonical_padloc(model_name: str, models: Path) -> List[str]:
    """
    Search in PADLOC models which models is canonical of another

    Args:
        model_name: name of the currently inspected model
        models: All other model in PADLOC model database

    Returns:
        List[str]: List of canonical models
    """
    canonical_sys = []
    if re.search("_other", model_name):
        basename = re.split("_other", model_name)[0]
        for canon_file in models.rglob("*.yaml"):
            name_canon = canon_file.stem
            if (re.search(f"{basename}_", name_canon) or re.search(f"{basename}$", name_canon)) \
                    and name_canon != model_name:
                canonical_sys.append(name_canon)
    return canonical_sys


def translate_padloc(padloc_db: Path, output: Path, binary_hmm: bool = False, hmm_coverage: float = None,
                     target_coverage: float = None, force: bool = False, disable_bar: bool = False) -> List[dict]:
    """
    Translate PADLOC models to PANORAMA models

    Args:
        padloc_db: Path to the PADLOC database
        output: Path to the directory where to write a HMM list file for PANORAMA annotation step
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        force: Flag to overwrite the output directory
        disable_bar: Flag to disable the progress bar

    Returns:
        List[dict]: List of dictionaries containing PANORAMA model information
    """
    list_data = []
    meta_df = parse_meta_padloc(padloc_db / "hmm_meta.txt")
    create_hmm_list_file(hmm_path=[padloc_db / "hmm"], output=output, metadata_df=meta_df, hmm_coverage=hmm_coverage,
                         target_coverage=target_coverage, binary_hmm=binary_hmm, force=force, disable_bar=disable_bar)
    logging.getLogger("PANORAMA").info("Begin to translate padloc models...")
    model_path = padloc_db / "sys"
    for model in tqdm(list(model_path.rglob("*.yaml")), unit="file", disable=disable_bar):
        canonical_sys = search_canonical_padloc(model.stem, model_path)
        data = read_yaml(model)
        list_data.append(translate_model_padloc(data, model.stem, meta_df, canonical_sys))
    return list_data


def translate_gene(elem: lxml.etree.Element, data: dict, hmm_df: pd.DataFrame, transitivity_mut: Callable) -> dict:
    """ Translate a gene from Defense-finder or MacSyFinder models into PANORAMA family models

    Args:
        elem: XML element corresponding to the gene to be translated
        data: Dictionary containing PANORAMA model information
        hmm_df: HMM information useful to translate gene into PANORAMA family models

    Returns:
        dict: Dictionary containing PANORAMA family model information
    """
    try:
        hmm_info = hmm_df.loc[elem.get('name')]
    except KeyError:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")
    else:
        if isinstance(hmm_info, pd.Series):
            dict_elem = {'name': hmm_info["protein_name"], 'presence': 'neutral', "parameters": {}}
            exchangeable_set = {hmm_info.name, hmm_info['secondary_name']}
        else:  # pd.Dataframe
            prot_name_list = hmm_info["protein_name"].mode().to_list()
            secondary_name_list = hmm_info["secondary_name"].mode().to_list()
            dict_elem = {'name': prot_name_list[0], 'presence': 'neutral', "parameters": {}}
            exchangeable_set = set(prot_name_list) | {elem.get('name')} | set(secondary_name_list)
        for attrib, value in elem.attrib.items():
            if attrib == 'presence':
                dict_elem['presence'] = value
            elif attrib == 'inter_gene_max_space':
                dict_elem['parameters']['transitivity'] = transitivity_mut(int(value))
            elif attrib == 'multi_system' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value == "True"):
                dict_elem["parameters"]["transitivity"] = -1
            elif attrib == 'multi_model' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_model'] = True

        for relation in elem:
            if relation.tag == "exchangeables":
                for gene in relation:
                    try:
                        exchangeable_info = hmm_df.loc[gene.get('name')]
                    except KeyError:
                        raise KeyError(f"In {data['name']}, one gene doesn't have a name")
                    else:
                        if isinstance(exchangeable_info, pd.Series):
                            exchangeable_set |= {exchangeable_info["protein_name"], exchangeable_info.name,
                                                 exchangeable_info['secondary_name']}
                        else:  # pd.Dataframe
                            prot_name_set = set(exchangeable_info["protein_name"].to_list())
                            secondary_name_set = set(exchangeable_info["secondary_name"].to_list())
                            name_set = set(exchangeable_info.index.to_list())
                            exchangeable_set |= prot_name_set | name_set | secondary_name_set
            else:
                logging.getLogger("PANORAMA").warning("Unexpected relation")
        if len(exchangeable_set) > 0:
            if isinstance(hmm_info, pd.Series):
                dict_elem["exchangeable"] = list(exchangeable_set.difference({hmm_info["protein_name"], ""}))
            else:  # pd.Dataframe
                prot_name_list = hmm_info["protein_name"].mode().to_list()
                dict_elem["exchangeable"] = list(exchangeable_set.difference({prot_name_list[0], ""}))
        return dict_elem


def translate_fu(elem: lxml.etree.Element, data: dict, hmm_df: pd.DataFrame, transitivity_mut: Callable):
    """ Translate a gene from Defense-finder or MacSyFinder models into PANORAMA functional unit models

    Args:
        elem: XML element corresponding to the gene to be translated
        data: Dictionary containing PANORAMA model information
        hmm_df: HMM information useful to translate gene into PANORAMA family models

    Returns:
        dict: Dictionary containing PANORAMA functional unit model information
    """
    try:
        name = elem.get('name')
    except KeyError:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")
    else:
        dict_elem = {'name': name, 'presence': 'neutral', "families": [],
                     'parameters': {}}

        for attrib, value in elem.attrib.items():
            if attrib == 'presence':
                dict_elem['presence'] = value
            elif attrib == 'min_mandatory_genes_required':
                dict_elem['parameters']['min_mandatory'] = int(value)
            elif attrib == 'min_genes_required':
                dict_elem['parameters']['min_total'] = int(value)
            elif attrib == 'inter_gene_max_space':
                dict_elem['parameters']['transitivity'] = transitivity_mut(int(value))
            elif attrib == 'multi_system' and value == 1:
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value is True):
                dict_elem["parameters"]["transitivity"] = -1
            elif attrib == 'multi_model' and value == 1:
                dict_elem['parameters']['multi_model'] = True

        for family in elem:
            dict_elem["families"].append(translate_gene(family, data, hmm_df))
        return dict_elem


def translate_macsyfinder_model(root: et.Element, model_name: str, hmm_df: pd.DataFrame,
                                canonical: List[str], transitivity_mut: Callable) -> dict:
    """
    Translate macsyfinder model into PANORAMA model

    Args:
        root: Model root information
        model_name: name of the model
        hmm_df: HMM information to translate models into PANORAMA models
        canonical: List of canonical models

    Returns:
        dict: PANORAMA model information
    """
    data_json = {"name": model_name, 'func_units': [], 'parameters': dict()}
    if len(canonical) > 0:
        data_json["canonical"] = canonical
    if root.attrib is not None:  # Read attributes root
        for parameter, value in root.attrib.items():
            if parameter == 'inter_gene_max_space':
                data_json['parameters']['transitivity'] = transitivity_mut(int(value))
            elif parameter == 'min_mandatory_genes_required':
                data_json['parameters']['min_mandatory'] = int(value)
            elif parameter == 'min_genes_required':
                data_json['parameters']['min_total'] = int(value)
            elif parameter == 'multi_system' and (value == 1 or value is True):
                data_json['parameters']['multi_system'] = True
            elif parameter == "loner" and (value == 1 or value is True):
                data_json["parameters"]["transitivity"] = -1
            elif parameter == 'multi_model' and (value == 1 or value is True):
                data_json['parameters']['multi_model'] = True

    fu_list = []
    fam_list = []
    for elem in root:
        if elem.tag == "gene":
            fam_list.append(translate_gene(elem, data_json, hmm_df, transitivity_mut))
        elif elem.tag == "functional_unit":
            fu_list.append(translate_fu(elem, data_json, hmm_df, transitivity_mut))
    if len(fu_list) == 0:  # only genes
        data_json["func_units"].append({'name': data_json["name"], 'presence': 'mandatory',
                                        "families": fam_list, 'parameters': data_json["parameters"]})
        data_json["parameters"] = {"transitivity": 0, "window": 1,  "min_mandatory": 1, "min_total": 1}
    else:
        for fu in fu_list:
            data_json["func_units"].append(fu)
    return data_json


def parse_macsyfinder_hmm(hmm: HMM, hmm_file: Path, panorama_acc: Set[str]) -> Dict[str, Union[str, int, float]]:
    """
    Read DefenseFinder HMM and get information for PANORAMA annotation step

    Args:
        hmm: HMM object
        hmm_file: Path to the DefenseFinder HMM
        panorama_acc: Set of PANORAMA accession ID for HMM

    Returns:
        List[Dict[str, Union[str, int, float]]]: List of dictionaries with all necessary information to write HMM metadata
    """
    hmm_dict = {"name": "", 'accession': "", "length": len(hmm.consensus), 'protein_name': "", 'secondary_name': "",
                "score_threshold": nan, "eval_threshold": 0.1, "ieval_threshold": 0.001, "hmm_cov_threshold": nan,
                "target_cov_threshold": nan, "description": ""}

    name = hmm.name.decode('UTF-8').split()
    if len(name) > 1:
        logging.getLogger("PANORAMA").error(f"HMM with problem to translate: {hmm_file.absolute().as_posix()}")
        raise IOError("It's not possible to get the name of the HMM. Please report an issue on GitHub.")
    else:
        if name[0] != hmm_file.stem.split('__')[-1]:
            hmm_dict["name"] = hmm_file.stem
            hmm_dict["protein_name"] = hmm_file.stem.split('__')[-1]
        else:
            hmm_dict["name"] = hmm_file.stem
            hmm_dict["protein_name"] = name[0]
    hmm.name = hmm_dict["protein_name"].encode('UTF-8')

    if hmm.accession is None or hmm.accession == "".encode("UTF-8"):
        hmm_dict["accession"] = gen_acc("PAN" + ''.join(choice(digits) for _ in range(6)), panorama_acc)
        hmm.accession = hmm_dict["accession"].encode("UTF-8")
    else:
        hmm_dict["accession"] = hmm.accession.decode("UTF-8")

    if hmm.description is not None:
        hmm_dict["description"] = hmm.description.decode("UTF-8")

    return hmm_dict


def create_macsyfinder_hmm_list(hmms_path: Path, output: Path, binary_hmm: bool = False, hmm_coverage: float = None,
                                target_coverage: float = None, force: bool = False,
                                disable_bar: bool = False) -> pd.DataFrame:
    """
    Read and parse all DefenseFinder HMM files and write a HMM list file for PANORAMA annotation step

    Args:
        hmms_path: Path to the HMM directory
        output: Path to the output directory where HMM list file will be written
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        force: Flag to overwrite the output directory
        disable_bar: Flag to disable the progress bar

    Returns:
        pd.Dataframe: Dataframe containing the HMM information needed for model translation
    """
    logging.getLogger("PANORAMA").info("Begin to translate HMM files from DefenseFinder...")
    hmm_dir = mkdir(output / 'hmm', force, erase=True)
    hmm_info_list = []
    panorama_acc = set()
    for hmm_file in tqdm(list(hmms_path.glob('*.hmm')), unit="HMM", desc='Translate HMM', disable=disable_bar):
        hmm_list = read_hmm(hmm_path=hmm_file)
        for hmm in hmm_list:
            hmm_dict = parse_macsyfinder_hmm(hmm, hmm_file, panorama_acc)
            hmm_dict["path"] = write_hmm(hmm, hmm_dir, binary_hmm).absolute().as_posix()
            hmm_info_list.append(hmm_dict)

    hmm_df = pd.DataFrame(hmm_info_list)
    hmm_df = hmm_df.sort_values(by=["name", "accession", "protein_name"], ascending=[True, True, True])
    hmm_df = hmm_df[['name', 'accession', 'path', 'length', 'protein_name', 'secondary_name', 'score_threshold',
                     'eval_threshold', "ieval_threshold", 'hmm_cov_threshold', 'target_cov_threshold', 'description']]
    if hmm_coverage is not None:
        logging.getLogger("PANORAMA").warning("HMM coverage threshold will be overwritten")
        hmm_df["hmm_cov_threshold"] = hmm_coverage
    if target_coverage is not None:
        logging.getLogger("PANORAMA").warning("target coverage threshold will be overwritten")
        hmm_df["target_cov_threshold"] = target_coverage
    hmm_df.to_csv(output / "hmm_list.tsv", sep="\t", index=False)
    hmm_df.set_index('name', inplace=True)
    logging.getLogger("PANORAMA").info("HMM list file created.")
    return hmm_df


def search_canonical_macsyfinder(model_name: str, models: Path) -> List[str]:
    """
    Search Canonical models for MacSyFinder models

    Args:
        model_name: Name of the model
        models: Path to all other models

    Returns:
        List[str]: List with name of canonical model
    """
    canonical_sys = []
    if re.search("-Type-", model_name):
        basename = re.split("-Type-", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        for canon_file in models.rglob(f"{base_class}*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}-Subtype-{base_type}-", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_sys.append(name_canon)
    elif re.search("_Type_", model_name):
        basename = re.split("_Type_", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        if len(base_type.split('_')) == 2:
            base_type, subtype = base_type.split('_')
        else:
            subtype = ""
        for canon_file in models.rglob("*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}_Type_{base_type}_{subtype}", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_sys.append(name_canon)
    elif model_name in ["CAS_Cluster", "CBASS", "Wadjet", "BREX"]:
        for canon_file in models.rglob(f"{model_name}*.xml"):
            if canon_file.stem != model_name:
                canonical_sys.append(canon_file.stem)
    return canonical_sys


def get_models_path(models: Path, source: str) -> Dict[str, Path]:
    """
    Associate a unique name for the model to its path in function of the source

    Args:
        models: Models database path
        source: Name of the source.

    Returns:
        Dict[str, Path]: A dictionary with model name for PANORAMA as key and the model path as value
    """
    model2path = {}

    for model in models.rglob("*.xml"):
        if source == "defense-finder":
            model2path[model.stem] = model
        elif source == "CONJScan":
            model2path[f"{model.parent.stem}_{model.stem}"] = model
        elif source == "TXSScan":
            model2path[model.stem] = model
        elif source == "TFFscan":
            model2path[model.stem] = model
        else:
            raise ValueError(f"Unknown source: {source}")

    return model2path


def translate_macsyfinder(macsy_db: Path, output: Path, binary_hmm: bool = False, hmm_coverage: float = None,
                          target_coverage: float = None, source: str = "", force: bool = False,
                          disable_bar: bool = False) -> List[dict]:
    """
    Translate MacSyFinder models into PANORAMA models and write all necessary file for PANORAMA steps

    Args:
        macsy_db: Models database path
        output: Path to output directory for PANORAMA files
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        source: Name of the source. Defaults to ""
        force: Flag to force overwrite files
        disable_bar: Flag to disable progress bar

    Returns:
         List[dict]: List of dictionaries containing translated models
    """
    def default_transitivity_mut(x: int):
        return x

    def CONJScan_transitivity_mut(x: int):
        return int(x/2) if int(x/2) >= x/2 else int(x/2) + 1

    if hmm_coverage is None:
        if source == "defense-finder":
            hmm_coverage = 0.4
        else:
            hmm_coverage = 0.5

    hmm_df = create_macsyfinder_hmm_list(macsy_db / "profiles", output, binary_hmm, hmm_coverage, target_coverage,
                                         force, disable_bar)
    list_data = []
    logging.getLogger('PANORAMA').info(f"Begin to translate {source} models")
    model_path = macsy_db / "definitions"
    if source == "CONJScan":
        transitivity_mut = CONJScan_transitivity_mut
    else:
        transitivity_mut = default_transitivity_mut
    for name, path in tqdm(get_models_path(model_path, source).items(), unit='file',
                           desc='Translate models', disable=disable_bar):
        try:
            canonical_sys = search_canonical_macsyfinder(name, model_path)
            root = read_xml(path)
            list_data.append(translate_macsyfinder_model(root, name, hmm_df, canonical_sys, transitivity_mut))
        except Exception as error:
            logging.getLogger('PANORAMA').error(f"Error translating {path.as_posix()}")
            raise Exception(f"Error translating {name} due to error: {error}")
    return list_data


def launch_translate(db: Path, source: str, output: Path,  binary_hmm: bool = False, hmm_coverage: float = None,
                     target_coverage: float = None, force: bool = False, disable_bar: bool = False):
    """
    Launch models translation process and write results for PANORAMA

    Args:
        db: Path to the models database that need to be translated
        source: Name of the source model. PADLOC, DefenseFinder and MacSyFinder are supported for now
        output: Path to the output directory to write all files needed for PANORAMA
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        force: Flag to force overwrite existing files
        disable_bar: Flag to disable progress bar
    """

    if source == "padloc":
        list_data = translate_padloc(padloc_db=db, output=output, binary_hmm=binary_hmm, hmm_coverage=hmm_coverage,
                                     target_coverage=target_coverage, force=force, disable_bar=disable_bar)
    elif source in ["defense-finder", "CONJScan", "TXSScan", "TFFscan"]:
        list_data = translate_macsyfinder(macsy_db=db, output=output, binary_hmm=binary_hmm, hmm_coverage=hmm_coverage,
                                          target_coverage=target_coverage, source=source,
                                          force=force, disable_bar=disable_bar)
    else:
        raise ValueError(f"The given source: {source} is not recognize. "
                         f"Please choose between padloc, defense-finder or macsy-finder")
    logging.getLogger("PANORAMA").info("Write models for PANORAMA...")
    model_list = []
    mkdir(output / 'models', force)
    for data in tqdm(list_data, unit="model", desc='Write models', disable=disable_bar):
        write_model(output, data)
        model_list.append([data['name'], output.resolve() / f"models/{data['name']}.json"])
    model_df = pd.DataFrame(model_list, columns=['name', 'path'])
    model_df = model_df.sort_values('name')
    if source == "CONJScan":
        plasmid_mask = model_df['name'].str.contains('Plasmids_')
        chromosome_mask = model_df['name'].str.contains('Chromosome_')
        model_df[plasmid_mask].to_csv(output / 'models_plasmids_list.tsv', sep="\t", header=False, index=False)
        model_df[chromosome_mask].to_csv(output / 'models_chromosome_list.tsv', sep="\t", header=False, index=False)
    model_df.to_csv(output / "models_list.tsv", sep="\t", header=False, index=False)
