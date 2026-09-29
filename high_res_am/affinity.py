"""
Top-level driver for an affinity maturation experiment.

AffinityMaturation owns the whole experiment: it turns the config into antigen/antibody
dictionaries, builds one GerminalCenter per `num_GCs`, and runs the injection series.

Call graph:
    run.py
      +-- AffinityMaturation.run()             (this file)
            +-- GerminalCenter.injection()     (lymph_node.py)  x num_GCs, per injection
                  +-- GerminalCenter.run_GC_cycle()             repeated until the GC ends

A GC finishes an injection when it goes extinct, exceeds its size cap, or depletes its
antigen pool below 20%. GCs that went extinct are dropped before the next injection, so
later injections run on fewer GCs than `num_GCs`.

Where to make changes:
  - How founders and antigens are assembled from the config -> build_Ab_Ag() below.
    Founder names must end in `_<integer>`: lymph_node.py parses that trailing integer as
    the parent_id used for per-founder ddG lookup in score.py.
  - Everything the GC actually does -> lymph_node.py.
"""

import os
import csv
import sys
import json
from timeit import default_timer as timer
from shutil import copyfile
from itertools import compress

# 10Jun24, HMP: importing functions and classes from other parts of this code repo
from lymph_node import BCell, GerminalCenter
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

        # 10Jun24, HMP: create logging directory and files; skip if debugging (debug==1)
        if self.config['experiment']['debug']==0:
            self.make_outputs() 
            self.save_config()

        self.dicAbs, self.dicAgs = self.build_Ab_Ag() # make dictionary of Abs and Ags
        self.make_scoring_function()
        self.Abs = list(self.dicAbs.keys()) 
        self.Ags = list(self.dicAgs.keys()) 
        self.n_injections = self.config['experiment']['num_injections']

        # Make the germinal centers (keep this last to avoid errors)
        self.gcs = self.make_gcs()
    
    def injection(self):
        dones = [gc.injection() for gc in self.gcs]
        return dones
    
    def run(self):
        """
        Run the affinity maturation algorithm
        """

        # Set up the simulation
        self.setup()
        
        start = timer()
        for i in range(self.n_injections):
            print("======= Starting injection {} of {} =======".format(i+1, self.n_injections))
            print("GCs: {}".format(len(self.gcs)))
            dones = self.injection()
            if all(dones):
                break
            survive = [False if i == True else True for i in dones]
            self.gcs = list(compress(self.gcs, survive))
            print("GCs that perished: {}".format(dones))
            print("GCs Post-Perish: {}".format(len(self.gcs)))
            print("======= Finished injection {} of {} =======".format(i+1, self.n_injections))
        end = timer()

        print('\nTotal time: %lfs, average per GC %lfs' % ((end - start), (end - start) / float(self.num_GCs)))
    

    # def finish(self):
    #     pass
    
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
        save the config.json file for reproducibility
        """

        if self.output_files is not None: # 10Jun24, HMP: this is ALWAYS true!
            print('Generating a copy of your configuration file, saved as config.json, in the /log directory for this experiment.') # hannah - revisit
            config_path = os.path.join(self.config['config_path'], 'config.json')
            with open(config_path, 'w') as f: 
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

        # JGF note - May need to change this for Hannah to do a concentration profile
        # Might not want to enumerate over the antigen_seqs, since they'll be the same.
        # Maybe could just take the index out of concentration
        for i, Ag in enumerate(self.config['antigen']['antigen_seqs']): # if not a cocktail, i = 0, Ag = the only sequence defined in "antigen_seqs"
            dicAgs = {
                self.config['antigen']['antigen_names'][i]:{'seq' : [Ag],
                                                            'conc' : self.config['antigen']['concentration']}
                                                            
                }
            
        # Build GL antibody
        dicAbs = dict()
        counts = dict()
        for Ab in self.config['antibody']['antibody_names']:
            if Ab not in counts:
                counts[Ab] = 0
            counts[Ab] += 1

            str_seq = ''
            len_of_each_region = dict()
            for region, seq in list(self.config['antibody']['antibody_seq'][Ab].items()):
                len_of_each_region[region] = len(seq)
                str_seq += seq

            dicAbs[Ab] = {
                    'seq' : str_seq,
                    'ratio' : None,
                    'lengths' : len_of_each_region,
                    'count' : counts[Ab]
                }
        

        return dicAbs, dicAgs