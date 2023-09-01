import os
import csv
import sys
from lymph_node import BCell, GerminalCenter
import json
from timeit import default_timer as timer
from score import scoring_function


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
        self.output_files = self.make_outputs()
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
        
        if self.config['log_path'] is None:
            return None
        elif not os.path.exists(self.config['log_path']):
            os.makedirs(self.config['log_path'])
        else:
            print("Log path already exists! Exiting...")
            sys.exit()
        
        #TODO: Make output files (need to check what we want to log)

    def make_gcs(self):
        gcs = [GerminalCenter(self.Abs,
                              self.Ags,
                              self.dicAgs,
                              self.dicAbs,
                              self.revert,
                              self.maturation_params,
                              self.score_config,
                              self.debug) for _ in range(self.num_GCs)]
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
        if self.output_files is not None:
            path = os.join.path(self.config['log_path'], 'config.json')

            with open(path, 'w') as f:
                json.dump(self.config, f, indent=3)


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
        for i, Ag in enumerate(self.config['antigen']['antigen_seqs']):
            dicAgs = {
                self.config['antigen']['antigen_names'][i]:{'seq' : [Ag],
                                                    'conc' : self.config['antigen']['concentration'][i]}
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
