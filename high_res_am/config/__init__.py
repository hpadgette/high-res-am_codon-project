import os
import json
from copy import deepcopy
import collections


def generic_config(pathogen):
    """
    Get the generic config
    """
    # with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "generic_config.json"), encoding='utf-8') as f:
    #     generic_config = json.load(f)

    # Get the default config for a given pathogen
    if pathogen in ["hiv"]:
        pathogen_config = "generic_hiv_no_comments.json"
    elif pathogen in ["influenza"]:
        pathogen_config = "generic_influenza.json"
    elif pathogen in ["malaria"]:
        pathogen_config = "generic_malaria.json"
    
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), pathogen_config), encoding='utf-8') as f:
        pathogen_config = json.load(f)
    generic_config = None
    print("generic_config: ", generic_config)
    return merge_config(generic_config, pathogen_config)

def exp_config(exp_config=None):
    """
    Get the experiment config
    """

    assert isinstance(exp_config, str), "Either no config provided, or incorrect formatting. Config must be a str."
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), exp_config), encoding='utf-8') as f:
        exp_config = json.load(f)

    try:
        # Check what pathogen was provided
        pathogen = exp_config["pathogen"]
        print("PATHOGEN: ", pathogen)
    except KeyError:
        # If no pathogen, default to hiv
        pathogen = "hiv"
        antibody = "VRC01"
        print("No pathogen provided in user config!\n \
               Using '{}' (BG505) with a single injection, and the '{}' GL as the default parameters.".format(pathogen, antibody))
    
    pathogen_config = generic_config(pathogen)
    pathogen_config = None
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
    print(generic_config)
    if generic_config is None:
        print("generic_config is None")
        exp_config['maturation_parameters']['p_CDR_affect'] = 1 - exp_config['maturation_parameters']['p_CDR_silent'] - exp_config['maturation_parameters']['p_CDR_lethal']
        exp_config['maturation_parameters']['p_exit'] = 1 - exp_config['maturation_parameters']['p_recycle']
        print(exp_config)
        return exp_config
    else:
        final_config = deepcopy(generic_config)
        for key, value in exp_config.items():

            # If the value isn't provided in the exp_config, get the default from generic_config
            if isinstance(value, collections.Mapping):
                final_config[key] = merge_config(generic_config.get(key, {}), value)
            else:
                final_config[key] = value

        final_config['maturation_parameters']['p_CDR_affect'] = 1 - final_config['maturation_parameters']['p_CDR_silent'] - final_config['maturation_parameters']['p_CDR_lethal']
        final_config['maturation_parameters']['p_exit'] = 1 - final_config['maturation_parameters']['p_recycle']
        return final_config