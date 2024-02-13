import os
import click
from config import get_exp_config
from datetime import datetime
from affinity import AffinityMaturation
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
    config = get_exp_config(config)
    
    if exp_name is not None:
        config['exp_name'] = exp_name

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
