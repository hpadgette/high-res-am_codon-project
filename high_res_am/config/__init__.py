"""
Configuration loading and merging.

A run's config is built in three layers, each overriding the one before it:

    config/generic_config.json          defaults shared by every run
    config/generic_<pathogen>.json      pathogen defaults, selected by the "pathogen" key
    your --config JSON                  the experiment itself

Two parameters are derived after merging and should not be set by hand:
    p_CDR_affect = 1 - p_CDR_silent - p_CDR_lethal
    p_exit       = 1 - p_recycle

JSON comment convention: json.load cannot parse `//` comments, so any file actually loaded
must stay comment-free. `generic_hiv.json` is the commented, human-readable reference;
`generic_hiv_no_comments.json` is the one this module reads. Edit one, edit both. Your own
experiment configs must likewise be plain JSON.

Where to make changes:
  - Adding a pathogen -> add generic_<pathogen>.json (and its no_comments twin) to this
    directory, then add the name to the branch in get_generic_config() below.
"""

import os
import json
from copy import deepcopy
import collections
import numpy as np

def get_generic_config(pathogen):
    """
    Get the generic config
    """
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "generic_config.json"), encoding='utf-8') as f:
        generic_config = json.load(f)

    # Get the default config for a given pathogen
    if pathogen in ["hiv"]:
        pathogen_config = "generic_hiv_no_comments.json"
    elif pathogen in ["influenza"]:
        pathogen_config = "generic_influenza.json"
    elif pathogen in ["malaria"]:
        pathogen_config = "generic_malaria.json"
    
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), pathogen_config), encoding='utf-8') as f:
        pathogen_config = json.load(f)
    return safe_merge_dicts(generic_config, pathogen_config)

def get_exp_config(exp_config=None):
    """
    Get the experiment config
    """

    assert isinstance(exp_config, str), "Either no config provided, or incorrect formatting. Config must be a str."
    with open(exp_config, encoding='utf-8') as f:
        exp_config = json.load(f)

    try:
        # Check what pathogen was provided
        pathogen = exp_config["pathogen"]
    except KeyError:
        # If no pathogen, default to hiv
        pathogen = "hiv"
        antibody = "VRC01"
        print("No pathogen provided in user config!\n \
               Using '{}' (BG505) with a single injection, and the '{}' GL as the default parameters.".format(pathogen, antibody))
    
    pathogen_config = get_generic_config(pathogen)
    
    final_config = safe_merge_dicts(pathogen_config, exp_config)
    final_config['maturation_parameters']['p_CDR_affect'] = np.round(1 - final_config['maturation_parameters']['p_CDR_silent'] - final_config['maturation_parameters']['p_CDR_lethal'], 
                                                                     decimals=3)
    final_config['maturation_parameters']['p_exit'] = np.round(1 - final_config['maturation_parameters']['p_recycle'], 
                                                               decimals=3)
    return final_config

def safe_merge_dicts(base_dict, update_dict):
    """
    Merges two dictionaries without changing the source dictionaries.
    Parameters
    base dict : dict
    Source dictionary with initial values.
    update dict : dict
    Dictionary with changed values to update the base dictionary.
    glatt
    Returns
    new_dict : dict
    Dictionary containing values from the merged dictionaries.
    """
    if base_dict is None:
        return update_dict
    base_dict = deepcopy(base_dict)
    for key, value in update_dict.items():
        if isinstance(value, collections.Mapping):
            base_dict[key] = safe_merge_dicts (base_dict.get(key, {}), value)
        else:
            base_dict[key] = value
    return base_dict