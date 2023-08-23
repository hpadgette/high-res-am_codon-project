import os
import json
from copy import deepcopy
import collections


def generic_config(pathogen):
    """
    Get the generic config
    """
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "generic_config.json"), encoding='utf-8') as f:
        generic_config = json.load(f)

    # Get the default config for a given pathogen
    if pathogen in ["hiv"]:
        pathogen_config = "generic_hiv.json"
    elif pathogen in ["influenza"]:
        pathogen_config = "generic_influenza.json"
    elif pathogen in ["malaria"]:
        pathogen_config = "generic_malaria.json"

    return merge_config(generic_config, pathogen_config)

def exp_config(exp_config=None):
    """
    Get the experiment config
    """

    assert isinstance(exp_config, str), "Either no config provided, or incorrect formatting. Config must be a str."
    exp_config = {}

    try:
        # Check what pathogen was provided
        pathogen = exp_config["pathogen"]
    except KeyError:
        # If no pathogen, default to hiv
        pathogen = "hiv (BG505_wt)"
        antibody = "VRC01"
        print("No pathogen provided in user config!\n \
               Using '{}' with a single injection, and the '{}' GL as the default parameters.".format(pathogen, antibody))
    
    pathogen_config = generic_config(pathogen)

    return merge_config(pathogen_config, exp_config)

def merge_config(generic_config, exp_config):
    """
    Merges the two config dictionaries

    Parameters
    ----------
        generic_config : dict
            generic config for a given pathogen

        exp_config : dict
            experimental setup for a given AM run
    
    Returns
    ----------
        final_config : dict
            configuration that will get fed to AM simulator
    """

    if generic_config is None:
        return exp_config
    final_config = deepcopy(generic_config)
    for key, value in exp_config.items():

        # If the value isn't provided in the exp_config, get the default from generic_config
        if isinstance(value, collections.Mapping):
            final_config[key] = merge_config(generic_config.get(key, {}), value)
        else:
            final_config[key] = value

    final_config['p_CDR_affect'] = 1 - final_config['p_CDR_silent'] - final_config['p_CDR_lethal']
    final_config['p_exit'] = 1 - final_config['p_recycle']
    return final_config