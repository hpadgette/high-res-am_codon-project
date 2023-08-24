import os
import click
from config import exp_config
from datetime import datetime
from mixture import AffinityMaturation
from copy import deepcopy

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
def main(config, exp_name):

    config = config if config != '' else None
    config = exp_config(config)
    
    if exp_name is not None:
        config['exp_name'] = exp_name

    timestamp = datetime.now().strftime("_%Y-%m-%d-%H%M%S")
    if config['experiment']['logging']:
        config['log_path'] = os.path.join('log/', config['exp_name'], timestamp)
    else:
        config['log_path'] = None

    run_AM(config)

