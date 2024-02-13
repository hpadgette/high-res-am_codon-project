import os
import csv
import sys
from lymph_node import BCell, GerminalCenter
import json
from timeit import default_timer as timer
from score import scoring_function
from shutil import copyfile


class AffinityMaturation:
    """
    Affinity Maturation Model environment

    Parameters
    ----------------------------------------
        config : dict
            configuration dictionary
    ----------------------------------------
    
    Attributes
    ----------------------------------------
    
    """

    def __init__(self, config=None):
        self.config = config
        self.set_config()

    def setup(self):
        
        # Create logging directory and files
        self.make_outputs()
        self.save_config()

        self.dicAbs, self.dicAgs = self.build_Ab_Ag()
        self.make_scoring_function()
        self.Abs = list(self.dicAbs.keys())
        self.Ags = list(self.dicAgs.keys())
        
        self.n_injections = self.config['experiment']['num_injections']
        
        # Make the germinal centers (keep this last to avoid errors)
        self.gcs = self.make_gcs()
    
    def injection(self):
        dones = [gc.injection() for gc in self.gcs]
        return all(dones)
    
    def run(self):
        """
        Run the affinity maturation algorithm
        """

        # Set up the simulation
        self.setup()
        
        start = timer()
        for i in range(self.n_injections):
            print("======= Starting injection {} of {} =======".format(i+1, self.n_injections))
            done = self.injection()
            if done:
                break
            print("======= Finished injection {} of {} =======".format(i+1, self.n_injections))
        end = timer()
        print('\nTotal time: %lfs, average per GC %lfs' % ((end - start), (end - start) / float(self.num_GCs)))
    

    def finish(self):
        pass
    
    def make_outputs(self):
        """
        Making the log path and setting the outputs
        """

        if self.config['experiment']["logging"]:
            for log in self.config['log_path']:
                if log is None:
                    continue
                elif not os.path.exists(log):
                    os.makedirs(log)
        else:
            print("Log path already exists! Exiting...")
            sys.exit()
        
        logs = ['config']
        for to_log in self.config['experiment']['log_level']:
            logs.append(to_log)
        
        self.output_files = logs
        

    def make_gcs(self):
        
        gcs = [GerminalCenter(self.Abs,
                              self.Ags,
                              self.dicAgs,
                              self.dicAbs,
                              self.revert,
                              self.maturation_params,
                              self.score_config,
                              self.config['log_path'],
                              i,
                              self.debug) for i in range(self.num_GCs)]
        return gcs
    
    def set_config(self):
        """
        Setting attributes from config to make it easier to read
        """

        self.num_GCs = self.config['experiment']['num_GCs']
        self.revert = self.config['experiment']['revert']
        self.score_config = self.config['experiment']['score_config']
        self.maturation_params = self.config['maturation_parameters']
        self.debug = self.config['experiment']['debug']

    def save_config(self):
        """
        Save the config.json file for reproducibility
        """

        # TODO: Need to change since we already set the log_path to a list

        if self.output_files is not None:
            config_path = os.path.join(self.config['config_path'], 'config.json')

            with open(config_path, 'w') as f:
                json.dump(self.config, f, indent=3)
        # if len(self.output_files) > 1:
        #     self.gc_log_paths = [os.path.join(self.config['log_path'], 'GC_{}_log.csv'.format(i)) 
        #                     for i in range(self.config['experiment']['num_GCs'])]
            


    def make_scoring_function(self):
        """
        Make the scoring function
        """
        self.scoring_function = scoring_function(self.score_config)
    
    def build_Ab_Ag(self):
        """
        Make the GL antibody and antigen
        """
        
        # Build antigens

        # May need to change this for Hannah to do a concentration profile
        # Might not want to enumerate over the antigen_seqs, since they'll be the same.
        # Maybe could just take the index out of concentration
        for i, Ag in enumerate(self.config['antigen']['antigen_seqs']):
            dicAgs = {
                self.config['antigen']['antigen_names'][i]:{'seq' : [Ag],
                                                    'conc' : self.config['antigen']['concentration']}
                }
        
        # Build GL antibody

        for i, Ab in enumerate(self.config['antibody']['antibody_names']):
            str_seq = ''
            len_of_each_region = dict()
            for region, seq in list(self.config['antibody']['antibody_seq'].items()):
                str_seq += seq
                len_of_each_region[region] = len(seq)
            

            dicAbs = {
                Ab : {
                    'seq' : str_seq,
                    'ratio' : None,
                    'lengths' : len_of_each_region
                }
            }
        return dicAbs, dicAgs
