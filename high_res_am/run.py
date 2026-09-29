"""
Entry point for an affinity maturation run.

    python -u -m run --config=config/my_experiment.json [--seed N] [--exp_name NAME]

Run this from the `high_res_am/` package directory. Mutability.csv and Substitution.csv are
opened by relative path during somatic hypermutation, so the working directory matters.

What this file does: load and merge the config, fix the RNG seed, build the log directory
tree, then hand off to AffinityMaturation.run() in affinity.py.

Where to make changes:
  - The system being studied (antigen, founder antibodies, ddG tables, model parameters)
    -> your experiment JSON in config/. No Python edits needed for any of this.
  - Which parameters exist and their defaults -> config/generic_config.json and
    config/generic_hiv_no_comments.json.
  - The GC reaction itself (selection, mutation, recycling) -> lymph_node.py.
  - How binding energy is computed -> score.py.

See the README for a reference table of every config parameter, including several that
appear in the example configs but currently have no effect.
"""

import os
import random
import click
import numpy as np
from config import get_exp_config
from datetime import datetime
from affinity import AffinityMaturation
from copy import deepcopy


def resolve_seed(config, cli_seed=None):
    """
    Decide the master RNG seed for this run, apply it, and record it in the config.

    Precedence: --seed on the command line > "seed" in the config > freshly drawn.

    Runs stay stochastic by default: with no seed given, one is drawn from OS entropy, so
    replicates differ from each other as intended. The seed that was actually used is
    written back into the config (and therefore into the saved config.json) and printed to
    stdout, so any individual run can be repeated later if it needs a closer look.
    """
    seed = cli_seed
    if seed is None:
        seed = config['experiment'].get('seed')
    if seed is None:
        seed = int(np.random.SeedSequence().generate_state(1)[0])

    seed = int(seed)
    np.random.seed(seed)
    random.seed(seed)
    config['experiment']['seed'] = seed
    return seed


def run_AM(config):
    """
    Run the affinity maturation algorithm
    """
    simulation = AffinityMaturation(deepcopy(config))
    
    # Run the simulation
    simulation.run()
    
    


@click.command()
@click.option('--config', default='', help='Path to the config file.')
@click.option('--exp_name', default=None, help='Name of the experiment for log path.')
@click.option('--seed', default=None, type=int,
              help='Master RNG seed. Omit to draw a fresh one (recorded in the log and the saved config).')
def main(config, exp_name, seed):

    config = config if config != '' else None
    config = get_exp_config(config)

    if exp_name is not None:
        config['exp_name'] = exp_name

    seed = resolve_seed(config, seed)
    print("======= RNG seed: {} =======".format(seed))
    print("To repeat this exact run: --seed {} (or add \"seed\": {} to the experiment block)".format(seed, seed))

    timestamp = datetime.now().strftime("_%Y-%m-%d-%H%M%S")
    if config['experiment']['logging']:
        config_path = os.path.join('log/', config['experiment']['exp_name'] + timestamp)
        config['config_path'] = config_path

        if len(config['experiment']['log_level']) > 0:
            # Make log paths for each specified log_level in the configuration
            if 'injection' in config['experiment']['log_level']:
                inj_log_path = os.path.join(config_path + '/injection_logs')
            else:
                inj_log_path = None

            if 'cycles' in config['experiment']['log_level']:
                cycle_log_path = os.path.join(config_path + '/cycle_logs')
            else:
                cycle_log_path = None

            if 'B Cell' in config['experiment']['log_level']:
                bCell_log_path = os.path.join(config_path + '/Bcell_logs')
            else:
                bCell_log_path = None
            
            log_path = [inj_log_path, cycle_log_path, bCell_log_path]
            config['log_path'] = log_path
            
        else:
            # If not just use the timestamp directory
            log_path = config_path
            config['log_path'] = log_path
    else:
        config['log_path'] = None


    run_AM(config)

if __name__ == '__main__': main()
