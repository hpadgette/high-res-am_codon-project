import os
import csv
import sys
from mixture import BCell, GerminalCenter
import json
from timeit import default_timer as timer

class AffinityMaturation():
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

    def setup(self):
        
        self.num_GCs = self.config['experiment']['num_gcs']
        # Create logging directory and files
        self.output_files = self.make_outputs()
        self.gcs = self.make_gcs()
        self.dicAbs, self.dicAgs = self.build_Ab_Ag()
        self.n_injections = self.config['experiment']['num_injections']
        self.debug = self.config['experiment']['debug']
        
    
    def injection(self):
        dones = [gc.injection(self.dicAgs, self.dicAbs) for gc in self.gcs]
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
        gcs = [GerminalCenter(body=self, debug=self.debug) for _ in range(self.num_GCs)]
        return gcs
    
    def set_config(self, config):
        pass

    def save_config(self):
        """
        Save the config.json file for reproducibility
        """
        if self.output_files is not None:
            path = os.join.path(self.config['log_path'], 'config.json')

            with open(path, 'w') as f:
                json.dump(self.config, f, indent=3)

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
            list_seq = list()
            len_of_each_region = dict()
            for region, seq in self.config['antibody']['antibody_seq'].values():
                list_seq.append(list(seq))
                len_of_each_region[region] = len(seq)
            str_seq = "".join(map(str, list_seq))

            dicAbs = {
                Ab : {
                    'seq' : str_seq,
                    'ratio' : None,
                    'lengths' : len_of_each_region
                }
            }
        return dicAbs, dicAgs
