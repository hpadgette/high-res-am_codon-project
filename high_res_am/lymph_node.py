import compute_descriptors_kayla
import csv
import sys
import modeller
import os
import importlib
import numpy as np  # numerical tools
from copy import deepcopy  # deepcopy copies a data structure without any implicit references
from timeit import default_timer as timer  # timer for performance
import importlib
import random
from Bio.Seq import Seq
import random
from collections import Counter
from collections import defaultdict
import warnings
from score import scoring_function
import h5py
warnings.filterwarnings("ignore")

np.random.seed()

class GerminalCenter(object):

    def __init__(self, Ab, Ags, dicAgs, dicAbs, 
                 revert, maturation_params,
                 score_config, log_paths, gc_id, debug=0):
        
        """ Initialize GC-specific variables. """
        self.debug = debug
        self.dicAgs = dicAgs
        self.dicAbs = dicAbs
        self.revert = revert
        self.maturation_params = maturation_params
        self.num_GC_founders = self.maturation_params['num_GC_founders']
        self.p_exit = self.maturation_params['p_exit']
        self.p_recycle = self.maturation_params['p_recycle']
        
        
        # Adding prob_metric for multiple epitopes
        self.prob_metric = self.maturation_params['prob_metric']
        
        self.Ab = Ab
        self.Ags = Ags
        self.energy_scale = self.maturation_params['energy_scale']
        self.n_injections = 0  # Number of injections so far
        self.cycle_number = 0
        self.nb_cycle_max = 99
        self.help_cutoff = self.maturation_params['help_cutoff']
        self.score_config = score_config
        self.parent_seq = self.dicAbs[Ab[0]]['seq']
        
        # Set fixed activation energy
        # self.activation_energy = self.maturation_params['activation_energy']
        
        # Set fixed delta between GL and E_a (hardcoded to be +0.75 kcal/mol for now)
        self.activation_energy = self.score_config['WT'] + 0.75


        self.gc_dict = deepcopy(self.__dict__)


        self.memory_cells = BCell(self.gc_dict, residues=None, nbs=1)
        self.memory_founders = BCell(self.gc_dict, residues=None, nbs=1)
        self.seed_nbs = self.maturation_params['pop_size']
        self.ag_start = self.maturation_params['ag_start']
        self.ag_left = self.ag_start
        self.p_scale = self.maturation_params['p_scale']
        residues = self.dicAbs[Ab[0]]['seq']
        self.B_cells = BCell(self.gc_dict, residues=residues, nbs=self.seed_nbs)
        self.gc_id = gc_id
        self.flag = self.maturation_params['flag']

        # Log paths should always be fed as a list regardless of the number of log levels
        print("ALL LOG PATHS: ", log_paths)
        for index, log_file in enumerate(log_paths):
            print(log_file)
            if "injection_logs" in log_file:
                print("MAKING INJECTION LOG")
                self.injection_log = open(log_file + '/GC_{}_log.csv'.format(self.gc_id), 'w')
                self.injection_log.write('injection,cycle_number,n_clones,nbs,avg_E,std_E,best_E,worst_E\n')
                self.injection_log.flush()
        
            elif "cycle_logs" in log_file:
                print("MAKING CYCLE LOG")
                self.cycle_log = open(log_file + '/GC_{}_log.csv'.format(self.gc_id), 'w')
                self.cycle_log.write('injection,cycle_number,n_clones,nbs,ag_left,avg_E,std_E,best_E,worst_E\n')
                self.cycle_log.flush()

            elif "Bcell_logs" in log_file:
                print("MAKING B CELL LOG")
                self.bcell_logPath = log_file + '/GC_{}_log.h5'.format(self.gc_id)
                self.bcell_log = h5py.File(log_file + '/GC_{}_log.h5'.format(self.gc_id), 'w')
                self.numeric_dataset = self.bcell_log.create_dataset('numeric_dataset', shape=(0, 6), maxshape=(None, 6), dtype='float64')
                self.numeric_dataset.attrs['description'] = "nbs, E, ag_captured, prob_capture, injection, cycle"
                self.sequence_dataset = self.bcell_log.create_dataset('sequence_dataset', shape=(0,), maxshape=(None,), dtype=h5py.string_dtype(encoding='ascii'))
            

        if isinstance(self.B_cells.nbs, int):
            self.GC_size_max = self.B_cells.nbs * 10
        else:
            self.GC_size_max = self.B_cells.nbs.sum() * 10

    def injection(self):
        self.n_injections += 1

        done = False

        # Set the antigens and concentrations
        self.Ag = self.dicAgs[self.Ags[0]]['seq']
        self.conc = self.dicAgs[self.Ags[0]]['conc']

        # Second injection onward, use memory founders
        if self.n_injections > 1:
            
            self.bcell_log = h5py.File(self.bcell_logPath, 'a')
            # Get numerical data, if it doesn't exist make a new dataset
            if 'numeric_dataset' in self.bcell_log:
                self.numeric_dataset = self.bcell_log['numeric_dataset']
            else:
                self.numeric_dataset = self.bcell_log.create_dataset('numeric_dataset', shape=(0, 6), maxshape=(None, 6), dtype='float64')
            
            # Get sequence data, if it doesn't exist make a new dataset
            if 'sequence_dataset' in self.bcell_log:
                self.sequence_dataset = self.bcell_log['sequence_dataset']
            else:
                self.sequence_dataset = self.bcell_log.create_dataset('sequence_dataset', shape=(0,), maxshape=(None,), dtype=h5py.string_dtype(encoding='ascii'))
        
            self.ag_left = self.ag_start 
            self.cycle_number = 0  # Reset cycle number
            self.plasma_cells = BCell(self.gc_dict, residues=None, nbs=0)  # Reset plasma cells

            # Set current B Cells to be memory founders
            self.B_cells = self.memory_founders
            print("Residues in Injection 2")
            print(self.B_cells.residues)
            # if len(self.B_cells.nbs) > 1 and not isinstance(self.B_cells.nbs, int):
                # print("Activation Energy for Injection {}: {}".format(self.n_injections, np.mean(self.B_cells.E)+1))
            # self.activation_energy = np.mean(self.B_cells.E) + 3.0
            # If there aren't any B Cells, exit since the GC has expired
            if isinstance(self.B_cells.nbs, int):
                if self.B_cells.nbs == 0:
                    if self.debug > 0:
                        print("GC failed; can't perform injection. Returning True")
                    return True
            elif len(self.B_cells.nbs) == 0:
                if self.debug > 0:
                    print("GC failed; can't perform injection. Returning True")
                return True

        # The value in this while loop should be a global var that is set at runtime
        while self.cycle_number < self.nb_cycle_max:
            self.run_GC_cycle()  # Run the GC cycle once
            """
            test_state() is checking to see if the current GC has exceeded the maximum
            size. It will also check if the maxmimum cycle_number has been reached. If
            neither flag is raised, it will then add 1 to the current cycle_number.
            """
            print("Done with GC Cycle!")
            if len(self.B_cells.nbs) > 1:
                print("Average deltaG for cycle {}: {}".format(self.cycle_number, np.average(self.B_cells.E, weights=self.B_cells.nbs)))
            print("There are {} Ags left for the next cycle!".format(self.ag_left))
            done, exceed_size_flag, ag_left_flag = self.test_state()
            print("Done? ", done)
            print("Too Big?? ", exceed_size_flag)
            print("Are we out of antigen? ", ag_left_flag)

            """
            JGF 1/30/24:
            I think we want to add the per-injection level logging here, and have it depend on any of the flags. There may 
            be an issue doing it this way because the GCs don't have any sort of identifier, so if we want logs for each GC
            we need to add this. We could also just have a column in the log file for the GC identifier? Either way I think
            we still need some sort of gc_id attribute.
            """
            # if not isinstance(self.B_cells.nbs, int):
            if len(self.B_cells.nbs) > 0:
                if self.injection_log is not None:
                    if done or exceed_size_flag or ag_left_flag:
                        print("We're writing the log!")
                        print("To be written\n")
                        print("{}, {}, {}, {}, {}, {}, {}, {} \n".format(self.n_injections, self.cycle_number, len(self.B_cells.nbs), self.B_cells.nbs.sum(), 
                            np.mean(self.B_cells.E), np.std(self.B_cells.E), np.min(self.B_cells.E), np.max(self.B_cells.E)))
                        self.injection_log.write('{},{},{},{},{},{},{},{}\n'.format(
                            self.n_injections, self.cycle_number, len(self.B_cells.nbs), self.B_cells.nbs.sum(), 
                            np.mean(self.B_cells.E), np.std(self.B_cells.E), np.min(self.B_cells.E), np.max(self.B_cells.E)))
                        self.injection_log.flush()
                        os.fsync(self.injection_log.fileno())
                if self.cycle_log is not None:
                    print("{}, {}, {}, {}, {}, {}, {}, {}, {} \n".format(self.n_injections, self.cycle_number, len(self.B_cells.nbs), 
                    self.ag_left, self.B_cells.nbs.sum(), np.mean(self.B_cells.E), np.std(self.B_cells.E), np.min(self.B_cells.E), np.max(self.B_cells.E)))
                    self.cycle_log.write('{},{},{},{},{},{},{},{},{}\n'.format(
                            self.n_injections, self.cycle_number, len(self.B_cells.nbs), self.B_cells.nbs.sum(), self.ag_left, 
                            np.mean(self.B_cells.E), np.std(self.B_cells.E), np.min(self.B_cells.E), np.max(self.B_cells.E)))
                    self.cycle_log.flush()
                    os.fsync(self.cycle_log.fileno())

            else:
                # When there aren't any B Cells left there's an error for the log writing; just writing np.inf when that happens
                if self.injection_log is not None:
                    self.injection_log.write('{},{},{},{},{},{},{},{}\n'.format(
                        self.n_injections, self.cycle_number, 0, np.inf, np.inf, np.inf, np.inf, np.inf))
                    self.injection_log.flush()
                    os.fsync(self.injection_log.fileno())
                if self.cycle_log is not None:
                    self.cycle_log.write('{},{},{},{},{},{},{},{}\n'.format(
                        self.n_injections, self.cycle_number, 0, np.inf, np.inf, np.inf, np.inf, np.inf))
                    self.cycle_log.flush()
                    os.fsync(self.cycle_log.fileno())

            if done:
                self.bcell_log.close()
                if self.debug > 1:
                    print("Finished after cycle {}".format(self.cycle_number))
                print('\n \n \n')
                break
            elif exceed_size_flag:
                self.bcell_log.close()
                self.memory_founders = self.pick_memCells()
                print('\n \n \n')
                break
            elif ag_left_flag:
                self.bcell_log.close()
                self.memory_founders = self.pick_memCells()
                print('\n \n \n')
                break
            print("Done with GC Cycle {}, on to the next one!".format(self.cycle_number))
            print('\n \n \n')

        return done

    def run_GC_cycle(self):
        print("=== Running GC Cycle ===")
        """ Run one cycle of the GC reaction. """

        self.run_dark_zone()  # DARK  ZONE - two rounds of division + SHM + updates cycle_number
        if self.debug > 2:
            print("After dark zone, {} B clones, {} total cells".format(len(self.B_cells.nbs), self.B_cells.nbs.sum()))
        total_cells = self.B_cells.nbs.sum()
        if total_cells <= 0 and self.debug > 0:
            print('GC extinct at cycle', self.cycle_number)
        else:
            self.run_binding_selection()  # LIGHT ZONE - selection for binding to Ag
            if self.debug > 2:
                print("After dark zone, {} B clones, {} total cells".format(len(self.B_cells.nbs),
                                                                            self.B_cells.nbs.sum()))
        total_cells = self.B_cells.nbs.sum()
        if total_cells <= 0:
            print('GC extinct at cycle', self.cycle_number)
        else:
            self.run_help_selection()  # LIGHT ZONE - selection to receive T cell help
            if self.debug > 2:
                print("After dark zone, {} B clones, {} total cells".format(len(self.B_cells.nbs),
                                                                            self.B_cells.nbs.sum()))
        total_cells = self.B_cells.nbs.sum()
        if total_cells <= 0:
            print('GC extinct at cycle', self.cycle_number)
        else:
            self.run_recycle()  # RECYCLE    - randomly pick exiting cells from the surviving B cells
            if self.debug > 2:
                print("After dark zone, {} B clones, {} total cells".format(len(self.B_cells.nbs),
                                                                            self.B_cells.nbs.sum()))
        if self.bcell_log is not None:
            print("Logging B Cells")
            timer_start = timer()                                                                           
            for i in range(len(self.B_cells.nbs)):
                floats_to_append = np.array([
                    self.B_cells.nbs[i], 
                    self.B_cells.E[i], 
                    self.B_cells.prob_capture[i], 
                    self.B_cells.ag_captured[i],
                    self.n_injections,
                    self.cycle_number])
                sequence_to_append = self.B_cells.residues[i][0]
                # print("RESIDUES FOR THE LOG: ", self.B_cells.residues[i][0])
                # print("floats: ", floats_to_append)
                # print("seq: ", sequence_to_append)
                self.append_data(floats_to_append, sequence_to_append)
            self.bcell_log.flush()
            timer_end = timer()
            
            print("Time to log B Cells: {}".format(timer_end - timer_start))
            
            
    def run_dark_zone_old(self, nb_rounds=2):
        """ B cells proliferate and undergo SHM in the dark zone. """
        for i in range(nb_rounds):
            self.B_cells.divide(self.cycle_number)
            self.B_cells.shm(self.Ag, self.cycle_number)
    
    def run_dark_zone(self, nb_rounds=2):
        print("=== Running Dark Zone ===")
        """ B cells proliferate and undergo SHM in the dark zone. """
        
        for i in range(nb_rounds):
            new_cells = []
            if isinstance(self.B_cells.nbs, int):
                self.calc_expansion_rate()
                self.B_cells.divide(self.cycle_number)
                self.B_cells.shm(self.Ag, self.cycle_number)
            else:
                all_ags_captured = [ag_captured for ag_captured in self.B_cells.ag_captured]
                self.calc_expansion_rate()
                self.B_cells.divide(self.cycle_number)
                self.B_cells.shm(self.Ag, self.cycle_number)

    def run_binding_selection_old(self):
        """ Select B cells for binding to antigen. """
        for i in range(len(self.B_cells.nbs)):
            Ag_bounds = np.exp(-1 * self.energy_scale * (self.B_cells.E[i] - self.activation_energy))  # Shape: (n_clones,)
            factors = self.conc * Ag_bounds
            langmuir_conjs = 1. / (1. + factors)
            n_dies = np.random.binomial(self.B_cells.nbs[i], langmuir_conjs)
            self.B_cells.nbs[i] -= n_dies
        self.B_cells.prune()
    
    def run_binding_selection(self):
        """ Select B cells for binding to antigen. """
        print("=== Binding Selection ===")
        new_cells=[]
        self.B_cells.prob_capture = np.zeros(shape=self.B_cells.nbs.shape)
        pvals = np.array([1./float(len(self.Ags))] * len(self.Ags))
        self.B_cells.last_bound = self.B_cells.vectorized_multinomial(ns=self.B_cells.nbs, condp = [1./float(len(self.Ags))] * len(self.Ags))[:]
        print("Ag Left at the start of binding selection: ", self.ag_left)
        print("ag_start at the start of binding selection: ", self.ag_start)
        for i in range(len(self.B_cells.nbs)):
            for j in range(len(self.Ags)):
                # compute binding energy and chance of death ( = 1 - chance of survival )
                Ag_bounds = np.exp(-1 * self.energy_scale * (self.B_cells.E[i] - self.activation_energy))  # Shape: (n_clones,)
                ag_conc       = self.ag_left/self.ag_start
                factor        = ag_conc * Ag_bounds
                langmuir_conj = 1. / (1. + factor)
                langmuir = 1 - langmuir_conj
                if langmuir - (ag_conc/(1+ag_conc))**ag_conc <= 0.01:
                    p_capture = 0
                else:
                    p_capture = (langmuir - (ag_conc/(1+ag_conc))**ag_conc) * self.p_scale
                self.B_cells.prob_capture[i] = p_capture
                
                # remove dead cells and update binding details
                n_die = np.random.binomial(self.B_cells.nbs[i], langmuir_conj)
                self.B_cells.nbs[i] -= n_die
                self.B_cells.last_bound[i] -= n_die
                
        self.B_cells.prune()
        self.B_cells.prob_capture = self.B_cells.prob_capture[self.B_cells.alive]
        self.assign_ag() # Explicit antigen capture
        
    
    def run_help_selection(self):
        print("=== Help Selection ===")
        """ Select B cells to receive T cell help. """
        binding_energies = self.B_cells.E
        
        # for i in range(len(self.B_cells.nbs)):
        #     binding_energies.append(self.B_cells.energy(self.Ag[0], self.cycle_number, self.B_cells.residues[i]))
        binding_energies = np.array(binding_energies)
        sort = np.argsort(binding_energies)[::-1]

        sorted_binding_energies = binding_energies[sort]
        sorted_nbs = self.B_cells.nbs[sort]
        
        nb_tot = self.B_cells.nbs.sum()
        if nb_tot > 0:
            cut_idx = int(self.help_cutoff * nb_tot) - 1
            energy_idx = np.argmax(np.cumsum(sorted_nbs) > cut_idx)
            # Returns index of first instance in which cumulative sum of nbs is greater than the cutoff value
            energy_cutoff = sorted_binding_energies[energy_idx]

            # Kill all B cells below threshold
            dies = binding_energies > energy_cutoff
            self.B_cells.nbs[dies] = 0

            # JGF 12/14/23: Something seems to be causing issues with the cut_idx because we don't prune, and
            # it's adding to the ties not taking away from them...
            
            # Pruning dead B Cells
            self.B_cells.prune()
            # Recomputing the cut_idx
            nb_tot = self.B_cells.nbs.sum()
            cut_idx = int(self.help_cutoff * nb_tot) - 1
            binding_energies = self.B_cells.E
            # Kill a fraction of B cells exactly at threshold
            # Tie deaths are split evenly across all clones exactly at threshold
            ties = binding_energies == energy_cutoff
            print("ties: ", ties)
            print("cut_idx: ", cut_idx)
            n_die_tie = self.B_cells.nbs.sum() - cut_idx
            print("n_die_tie: ", n_die_tie)
            # Total number of cells to remove for clones exactly at threshold
            fraction_die_tie = n_die_tie / self.B_cells.nbs.sum()
            print("fraction_die_tie: ", fraction_die_tie)
            print("total nbs pre tiedie: ", self.B_cells.nbs)
            print("nbs pre tiedie: ", self.B_cells.nbs[ties])
            # Fraction of cells to remove for each clone exactly at threshold
            self.B_cells.nbs[ties] = np.multiply(self.B_cells.nbs[ties], (1 - fraction_die_tie), casting="unsafe")
            print("nbs post tiedie: ", self.B_cells.nbs[ties])
            # Note: Can't use *= because we're multiplying an int by a float
            print("multiplication for tiedie", np.multiply(self.B_cells.nbs[ties], (1 - fraction_die_tie), casting="unsafe"))
            self.B_cells.prune()

    def run_recycle(self):
        print("=== Recycling ===")
        """ Randomly select B cells to be recycled back into the GC or to become memory cells. """
        
        n_exits = np.random.binomial(self.B_cells.nbs, self.p_exit)
        self.B_cells.nbs -= n_exits
        new_memory_cells = BCell(self.gc_dict, residues=self.B_cells.residues, nbs=n_exits)
        new_memory_cells.E = self.B_cells.E
        self.B_cells.prune()

        if self.revert:  # Old behavior: Don't add memory cells for the last cycle
            GC_size = self.B_cells.nbs.sum()
            if (self.cycle_number == self.nb_cycle_max - 1) or (GC_size > self.GC_size_max):
                pass
            else:
                self.memory_cells.clone(new_memory_cells.residues, new_memory_cells.nbs)  # Always add memory cells
        else:
            # TODO: (JGF) I think this should be replaced with appending the memory cells each time this
            # function gets called. Otherwise, the memory cells only can be sampled from the last iteration.
            # I would think that it's more accurate to have most of the memory cells pulled from the point in AM
            # where they look similar to the germline. Plus, then we don't have to tune how many GL vs mem_cells.
            # Need to add an append_mem_cells() function I think, or add it to this.
            if new_memory_cells.nbs.sum() > 0:
                self.memory_cells.append_mem_cells(new_memory_cells)  # Always add memory cells

    def pick_memCells(self):
        print("=== Picking Memory Cells ===")
        # For simplicity, memory cells are chosen with replacement.
        # To do it without replacement, first choose from self.memory_cells.nbs.sum(),
        # then find the clones corresponding to each chosen ID. Or, loop what's below
        # selecting 1 cell at a time, then reduce the corresponding nbs by 1 each time.
        nbs = self.memory_cells.nbs
        if isinstance(nbs, int):
            n_mem_clones = nbs
            nb_tot = nbs
        else:
            n_mem_clones = len(nbs)  # Number of unique clones
            nb_tot = nbs.sum()

        if self.flag == 1:
            i = 0
        else:
            i = 2

        if self.revert:
            if n_mem_clones < self.num_GC_founders:
                if self.debug > 0:
                    print("Not enough memory cells to form founders!")
                return BCell(self.gc_dict, residues=None, nbs=0)
            p = None
            replace = False
        else:
            if nb_tot < self.num_GC_founders:
                if self.debug > 0:
                    print("Not enough memory cells to form founders!")
                return BCell(self.gc_dict, residues=None, nbs=0)
            p = nbs / nb_tot # Probabilities for choosing each clone (proportaional to nbs)
            replace = True # Approximation
        print("Number of Memory Clones: ", n_mem_clones)
        id_new_founders = np.random.choice(n_mem_clones, self.num_GC_founders, p=p, replace=replace)
        print("ID of New Founders: ", id_new_founders)

        #TODO: Allow for multiple GC founders after first injection
        residues = self.memory_cells.residues[id_new_founders][:]
        
        if self.revert:
            nbs = self.memory_cells.nbs[id_new_founders]
        else:
            # nbs = self.seed_nbs  # Instantly proliferate to 128
            nbs = 10 * np.ones(self.num_GC_founders, dtype=np.int) # JGF: Testing multiple memory cells
        
        new_founders = BCell(self.gc_dict, residues=residues, nbs=nbs)
        new_founders.E = self.memory_cells.E[id_new_founders]
        return new_founders

    def test_state(self):
        print("===Testing State===")
        GC_size = self.B_cells.nbs.sum()  # total number of cells in the GC
        print("Size of GC: ", GC_size)

        if (self.cycle_number == self.nb_cycle_max - 1) or (
                GC_size > self.GC_size_max):  # at the end, all B cells exit the GC
            self.plasma_cells = self.B_cells

        if (GC_size <= 0):
            exceed_size_flag = False
            ag_left_flag = False
            done = True
            if self.debug > 0:
                print("Done: GC size = 0 (all dead)")

        elif (GC_size > self.GC_size_max):
            exceed_size_flag = True
            done = False
            ag_left_flag = False

        elif self.ag_left/self.ag_start < 0.2:
            ag_left_flag = True
            done = False
            exceed_size_flag = False

        else:
            exceed_size_flag = False
            done = False
            ag_left_flag = False

        self.cycle_number += 1
        
        return done, exceed_size_flag, ag_left_flag
    
    def assign_ag(self):
        print("===Assigning Ag===")
        """Assign Ag to B cell clones based on Langmuir probability
        Each B cell is assumed to have 100 BCRs which can each extract 
        up to 2 Ag molecules yielding a total of 200 Ag molecules maximum 
        per B cell. (Amitai et al. PLOS Comp. Bio. 2018)"""
        
        print("\n")
        assert len(self.B_cells.prob_capture) == len(self.B_cells.nbs)
        # live_clones = [(p,b) for p, b in zip(self.prob_capture, self.B_cells) if self.B_cells.nbs > 0]
        # try:
        #     p_sorted = sorted(self.B_cells.prob_capture, reverse=True)
        # except:
        #     return self.ag_left

        # JGF: Was enumerate(p_sorted)
        for i in np.argsort(-1*self.B_cells.prob_capture):
            p = self.B_cells.prob_capture[i]
            
            """
            JGF: 
            Ag to capture is dependent on the number of antigens present in the GC, but
            should this be based on the max number each clone can get? If we use larger values then each
            clone wants thousands of Ag copies, but there would probably be way more than 1000 Ags in the GC?
            """
            # chances = int(np.round(self.ag_left/200, 0))
            chances = 200
            # print("Chances for clone {}: {}".format(i, chances))
            ag_to_capture = np.random.binomial(chances, p)
            wanted_by_clone = self.B_cells.nbs[i]*ag_to_capture
            print("Energy: ", self.B_cells.E[i])
            print("prob_capture: ", np.round(p, 4))
            print("Ag clone wants to capture: ", ag_to_capture)
            print("Ag already captured: ", self.B_cells.ag_captured[i]," per B cell.")
            print("Clone has ", self.B_cells.nbs[i], " B cells.")
            print("Total wanted by clone: ", wanted_by_clone)
            if wanted_by_clone < self.ag_left:
                if self.B_cells.ag_captured[i] + ag_to_capture <= 200:
                    self.B_cells.ag_captured[i] += ag_to_capture
                    self.ag_left -= wanted_by_clone
                    print("Assigned", ag_to_capture, "to each B cell of clone", i, ". There is", self.ag_left, "left.")
                else:
                    ag_needed = 200 - self.B_cells.ag_captured[i]
                    self.B_cells.ag_captured[i] = 200
                    self.ag_left -= ag_needed
                    print("B cell reached max ag!")
            elif self.ag_left > 0:
                divided_amount = np.floor(self.ag_left/self.B_cells.nbs[i])
                if self.B_cells.ag_captured[i] + divided_amount <= 200:
                    print("Not enough Ag left. Assigning remaining", self.ag_left, "as", divided_amount, "to each B cell.")
                    ag_needed = divided_amount * self.B_cells.nbs[i]
                    self.B_cells.ag_captured[i] += divided_amount
                    self.ag_left -= ag_needed
                    print("Assigned", divided_amount, "to each B cell of clone", i, ". There is", self.ag_left, "left.")
                else:
                    ag_needed = 200 - self.B_cells.ag_captured[i]
                    self.B_cells.ag_captured[i] = 200
                    self.ag_left -= ag_needed
                    print("B cell reached max ag of", self.B_cells.ag_captured[i], "!")
            else:
                print("No Ag left!")
            print("\n")
    
    def calc_expansion_rate(self):
        print("===Calculating Expansion Rate===")
        """
        JGF (12/8/23):
        Testing a different implementation. The current one basically gaurantees the first to 30 Ags wins, and
        gives it a div_rate of 3. Additionally, if a clone as 200 or 30, they still get the same expansion bonus.
        Experiments show max should be 6, and I would assume it would follow some sort of Michaelis-Menten sort of pattern.
        
        Testing with: div_max*(n_cap/n_max)/(1/(A*n_cap+B) + n_cap/n_max) # A and B are empirical & subject to change
        """

        div_max = 6
        n_max = 200

        for i, n_cap in enumerate(self.B_cells.ag_captured):
            if n_cap > 0:
                self.B_cells.div[i] = np.round((div_max*(n_cap/n_max))/((1/(0.1*n_cap + 0.1)) + (n_cap/n_max)) + 0.5, 0)
            else:
                self.B_cells.div[i] = 1
            
            print("Antigens Captured: {} \n Division Rate: {}".format(n_cap, self.B_cells.div[i]))

        # for i, ag in enumerate(self.B_cells.ag_captured):
        #     if ag <= 30:
        #         self.B_cells.div[i] = 1
        #     else:
        #         if ag >= np.percentile(self.B_cells.ag_captured, 95):
        #             self.B_cells.div[i] = 3
        #         elif ag >= np.percentile(self.B_cells.ag_captured, 80):
        #             self.B_cells.div[i] = 2
        #         else:
        #             self.B_cells.div[i] = 1
    
    def append_data(self, numeric_data, sequence_data):
        # Append numeric data
        current_numeric_length = self.numeric_dataset.shape[0]
        new_numeric_length = current_numeric_length + 1
        self.numeric_dataset.resize((new_numeric_length, 6))
        self.numeric_dataset[current_numeric_length] = numeric_data

        # Append string data
        current_string_length = self.sequence_dataset.shape[0]
        new_string_length = current_string_length + 1
        self.sequence_dataset.resize((new_string_length,))
        self.sequence_dataset[current_string_length] = sequence_data.encode('ascii')

class BCell():
    
    def __init__(self, gc_dict, residues, nbs, **kwargs):
        """ Initialize clone-specific variables. 
            nbs          - population size
            residues         - array of the values for each residue of the B cell
            E           - total binding energy
            nb_CDR_mut  - number of accumulated CDR mutations
            antigens    - list of antigens
            breadth     - breadth against test panel
            nb_Ag       - number of antigens available
            last_bound  - number of individuals that last bound each Ag
            mut_res_id  - index of mutated residue
            delta_res   - incremental change at residue
            generation  - generation in the GC reaction
            cycle_number - cycle number
            history     - history of mutations, generation occurred, and effect on Q/Ec and breadth against test panel """
        # NOTE: this is the number of mutations in the GL CODON sequence not AA sequence (# aa is < or = # nt mutations)
        if 'memory_cells' in gc_dict:
            gc_dict.pop('memory_cells')
        elif 'memory_founders' in gc_dict:
            gc_dict.pop('memory_founders')
        #TODO: I think this is bugged, doesn't seem to be loading in the gc_dict as expected
        for k, v in gc_dict.items():
            self.__dict__[k] = v
        self.residues = residues
        
        self.nbs = nbs
        if 'maturation_params' in gc_dict:
            self.p_mut = self.maturation_params['p_mut']
            self.p_CDR_affect = self.maturation_params['p_CDR_affect']
            self.p_CDR_silent = self.maturation_params['p_CDR_silent']
            self.p_CDR_lethal = self.maturation_params['p_CDR_lethal']
            self.p_CDR = self.maturation_params['p_CDR']
            self.p_metric = self.maturation_params['prob_metric']

        # Conditional p values used in vectorized_multinomial() are repeated often, so compute them once here.
        p = np.array([self.p_CDR_lethal, self.p_CDR_silent, self.p_CDR_affect])
        ps = p.cumsum(axis=-1)
        with np.errstate(divide='ignore', invalid='ignore'):
            condp = p / ps
        condp[np.isnan(condp)] = 0.0
        self.shm_condp = condp
        self.scoring_function = scoring_function(self.score_config).__call__
        
        if 'nb_CDR_mut' in kwargs:
            self.nb_CDR_mut = kwargs['nb_CDR_mut']
        else:
            self.nb_CDR_mut = 0

        if 'cycle_number' in kwargs:
            self.cycle_number = kwargs['cycle_number']
        else:
            self.cycle_number = 2

        if 'E' in kwargs:
            self.E = kwargs['E']
        else:
            if isinstance(nbs, int) and isinstance(residues, str):
                self.E = self.energy(self.dicAgs[self.Ags[0]]['seq'][0], self.cycle_number, self.residues)
            # elif residues is not None:
            #     E = np.zeros(shape=(len(residues),))
            #     for i, res in enumerate(residues):
            #         E[i] = self.energy(self.dicAgs[self.Ags[0]]['seq'][0], self.cycle_number, res)
            #     self.E = E
            else:
                self.E = np.inf

        if 'generation' in kwargs:
            self.generation = kwargs['generation']
        else:
            self.generation = 0
        
        self.E = np.array(self.E).reshape((1,))
        if isinstance(self.nbs, int):
            self.ag_captured = 0
            self.ag_captured = np.array(self.ag_captured).reshape((1,))
            self.div = np.array([1]).reshape((1,))
        else:
            print("NBS SHAPE FOR NEW B CELLS: ", self.nbs.shape)
            print("RESIDUES SHAPE FOR NEW B CELLS: ", self.residues.shape)
            self.ag_captured = np.zeros(self.nbs.shape[0], dtype=int)
            self.div = np.ones(self.nbs.shape[0], dtype=int)
        
        if gc_dict['maturation_params']['split_founders'] and isinstance(self.nbs, int):
            nbs = int(np.round(nbs/self.maturation_params['num_GC_founders']))
            self.nbs = np.repeat(nbs, self.maturation_params['num_GC_founders'])
            # print("NBS INIT: ", self.nbs)
            # print("LENGTH NBS: ", len(self.nbs))
            self.E = np.repeat(self.E[0], self.maturation_params['num_GC_founders'])
            # print("E INIT: ", self.E)
            # print("LENGTH E: ", len(self.E))
            self.ag_captured = np.repeat(0, self.maturation_params['num_GC_founders'])
            # print("AG INIT: ", self.ag_captured)
            # print("LENGTH AG: ", len(self.ag_captured))
            self.div = np.repeat(1, self.maturation_params['num_GC_founders'])
            # print("DIV INIT: ", self.div)
            # print("LENGTH DIV: ", len(self.div))
            for i in range(0, self.maturation_params['num_GC_founders']-1):
                self.residues = np.vstack([self.residues, residues])
            # print("LENGTH RESIDUES: ", len(self.residues))
        # self.prob_capture = np.array(0).reshape((1,))
        self.lenFRWH1 = self.dicAbs[self.Ab[0]]['lengths']['FRWH1_nt']
        self.lenCDRH1 = self.dicAbs[self.Ab[0]]['lengths']['CDRH1_nt']
        self.lenFRWH2 = self.dicAbs[self.Ab[0]]['lengths']['FRWH2_nt']
        self.lenCDRH2 = self.dicAbs[self.Ab[0]]['lengths']['CDRH2_nt']
        self.lenFRWH3 = self.dicAbs[self.Ab[0]]['lengths']['FRWH3_nt']
        self.lenCDRH3 = self.dicAbs[self.Ab[0]]['lengths']['CDRH3_nt']
        self.lenFRWL1 = self.dicAbs[self.Ab[0]]['lengths']['FRWL1_nt']
        self.lenCDRL1 = self.dicAbs[self.Ab[0]]['lengths']['CDRL1_nt']
        self.lenFRWL2 = self.dicAbs[self.Ab[0]]['lengths']['FRWL2_nt']
        self.lenCDRL2 = self.dicAbs[self.Ab[0]]['lengths']['CDRL2_nt']
        self.lenFRWL3 = self.dicAbs[self.Ab[0]]['lengths']['FRWL3_nt']
        self.lenCDRL3 = self.dicAbs[self.Ab[0]]['lengths']['CDRL3_nt']
        # if self.residues is not None and isinstance(self.residues, str):
        # print("Parent Sequence: ", self.parent_seq)
        # self.parent_seq = deepcopy(str(Seq(self.parent_seq).translate()))
        # self.parent_seq = self.parent_seq[0:118] + '/' + self.parent_seq[118:]
        self.len_heavy_nt = self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3
        self.len_light_nt = self.lenFRWL1 + self.lenCDRL1 + self.lenFRWL2 + self.lenCDRL2 + self.lenFRWL3 + self.lenCDRL3
        self.len_heavy_aa = len(str(Seq(self.parent_seq[0:self.len_heavy_nt]).translate()))
        self.len_light_aa = len(str(Seq(self.parent_seq[self.len_heavy_nt+1:]).translate()))
        self.parent_seq = deepcopy(str(Seq(self.parent_seq).translate()))
        self.parent_seq = self.parent_seq[0:len_heavy_aa] + '/' + self.parent_seq[len_heavy_aa+1:]


    def clone(self, residues, nbs):
        print("===Cloning B Cells===")
        """ Add clones with given residues and nbs """
        
        # If nbs is an int, repeat that value for all new clones
        if isinstance(nbs, int):
            nbs = np.repeat(nbs, residues.shape[0])
            E = np.repeat(self.E[0], residues.shape[0])
            ag_captured = np.repeat(0, residues.shape[0])
            div = np.repeat(1, residues.shape[0])
        # print("nbs in Clone(): ", len(self.nbs))
        # print("nbs to add in Clone(): ", len(nbs))
        # print("Num AG Captured in Clone(): ", len(self.ag_captured))
        # print("Num Ag Captured to add in Clone(): ", len(ag_captured))
        self.residues = np.vstack([self.residues, residues])
        self.nbs = np.concatenate([self.nbs, nbs])
        self.E = np.concatenate([self.E, E])
        self.ag_captured = np.concatenate([self.ag_captured, ag_captured])
        self.div = np.concatenate([self.div, div])


    # Might be unnecessary
    def update_history(self):
        """ Add current parameters to the history list. """
        self.history['generation'].append(self.generation)
        self.history['residues'].append(self.residues)
        self.history['nb_CDR_mut'].append(self.nb_CDR_mut)
        self.history['mut_res_id'].append(self.mut_res_id)
        self.history['E'].append(self.E)
        self.history['delta_res'].append(self.delta_res)
        self.history['cycle_number'].append(self.cycle_number)


    def energy(self, Ag, cycle_number, residues):
        
        """ Return binding energy with input antigen. """
        
        # Choosing epitope to bind based on list of probabilities in p_metric
        # This should return the index of the choice e.g. for p_metric = [0.9, 0.1] there is a 90% chance
        # it will return 0, and a 10% chance it returns 1
        if isinstance(self.p_metric, list):
            binding_site = random.choices(range(len(self.p_metric)), weights=self.p_metric, k=1)[0]
            print("Binding Site Index: ", binding_site)
        else:
            binding_site = 0
        
        # Convert to Bio nt_seq
        Ab_nt = Seq(residues)

        # Convert to AA
        Ab_aa = str(Ab_nt.translate())

        # Format for chains in Modeller
        Ab = Ab_aa[0:118] + '/' + Ab_aa[118:]
        if "*" not in Ab_aa:
            if self.score_config['scoring_function'] == 'singlepoint':
                dg = self.scoring_function(ab=Ab, wt_ab=self.parent_seq)
            else:
                print("Ag/Ab seq to scoring funcion: {}/{}".format(Ag, Ab))
                dg = self.scoring_function(self.scoring_function, ag=Ag, ab=Ab, binding_site=binding_site)
        else:
            dg = 100
        # print("deltaG: ", dg)
        return dg

    
    def divide(self, cycle_number):
        print("=== dividing ===")
        """ Run one round of division. """   
        #BMP change to account for differential expansion rates (Gitlin et al. 2014)

        # if cycle_number == 2:
        #     print("Dividing 1 time in round", cycle_number)
        #     self.nbs *= 2
        #     self.generation += 1
        #     self.cycle_number = cycle_number
        # else:
        if len(self.div) < 1:
            print("Everyone is dead, we can't divide!")
            print("======= Exiting =======")
            sys.exit()

        print("Dividing", self.div, "times in round", cycle_number)
        pop_growth = [int(2*i) for i in self.div]
        new_ag_captured = [int(ag_captured/(2*i)) for ag_captured,i in zip(self.ag_captured,self.div)]
        self.nbs = np.multiply(self.nbs, pop_growth)
        self.ag_captured = new_ag_captured
        self.generation += 1
        self.cycle_number = cycle_number


    def calculate_breadth(self, testpanel, threshold, panelSize):
        test_energies = np.inner(testpanel, self.residues)
        return np.mean(test_energies > threshold, axis=0)


    def mutate_CDR(self, Ag, n_clones):
        print("=== Mutating CDR ===")

        """
        Change in energy due to affinity-affecting CDR mutation. Only one residue mutates.

        Targeting model: choose nt to mutate based on 5mer mutability and sequence occurrences

        Substitution model: choose nt to mutate to based on probability of middle nt in 5mer
                            mutating to any other nt
        """
        start1 = timer()
        global picked_list2
        Fivemer = []
        Fivemer_index = []
        
        residues = self.residues[-n_clones:]
        cln = 0
        for clone_residues in residues:
            clone_residues = clone_residues[0]
            for i in range(len(str(clone_residues))):  # Iterating over the entire sequence
                if i < (len(str(clone_residues)) - 4):
                    # If i is less than four spaces from the end of the sequence
                    # (i.e., there is a fivemer available)
                    if self.lenFRWH1 <= (i + 2) < self.lenFRWH1 + self.lenCDRH1:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                        # Appending the section of the sequence to check for mutability

                    elif self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 <= (i + 2) < self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 <= (
                            i + 2) < self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 <= (
                            i + 2) < self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 + self.lenCDRL1:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 + self.lenCDRL1 + self.lenFRWL2 <= (
                            i + 2) < self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 + self.lenCDRL1 + self.lenFRWL2 + self.lenCDRL2:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 + self.lenCDRL1 + self.lenFRWL2 + self.lenCDRL2 + self.lenFRWL3 <= (
                            i + 2) < self.lenFRWH1 + self.lenCDRH1 + self.lenFRWH2 + self.lenCDRH2 + self.lenFRWH3 + self.lenCDRH3 + self.lenFRWL1 + self.lenCDRL1 + self.lenFRWL2 + self.lenCDRL2 + self.lenFRWL3 + self.lenCDRL3:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])

            reader = csv.DictReader(open('Mutability.csv'))
            mutTable = [obs for obs in reader]
            a = [float(o['Mutability']) for o in mutTable]
            mutSum = np.sum(a)

            probs = []
            for i in Fivemer:
                FivemerMutTable = self.selectObservations(mutTable, ['Fivemer'], [str(i)])
                a = [float(o['Mutability']) for o in FivemerMutTable]
                probs.append(a[0])
            # print("probs: ", probs)
            probsNorm = probs
            probsNorm /= np.sum(probsNorm)
            
            picked_int = np.random.choice(range(len(Fivemer)), 1, p=probsNorm)
            picked_5mer = Fivemer[picked_int[0]]
            picked_5mer_index = Fivemer_index[picked_int[0]]

            reader2 = csv.DictReader(open('Substitution.csv'))
            mutTable2 = [obs for obs in reader2]
            bases = ['A', 'C', 'G', 'T']
            probs2 = []
            temp = []
            temp.append(picked_5mer)

            for i in temp:
                FivemerMutTable2 = self.selectObservations(mutTable2, ['Fivemer'], [str(i)])
                A = [float(o['A']) for o in FivemerMutTable2]
                C = [float(o['C']) for o in FivemerMutTable2]
                G = [float(o['G']) for o in FivemerMutTable2]
                T = [float(o['T']) for o in FivemerMutTable2]
                probs2.extend([A, C, G, T])
                probs3 = np.concatenate(probs2, axis=0)
                picked_list2 = np.random.choice(bases, 1, p=probs3)
            picked_ntmut = picked_list2[0]

            temp_res1 = list(deepcopy(str(clone_residues)))
            index = picked_5mer_index
            self.mut_res_id = index
            delta = picked_ntmut
            self.delta_res = delta
            # print("index/delta:", index, delta)
            
            temp_res1 = [i for i in temp_res1 if i == 'A' or i == 'T' or i == 'G' or i == 'C']
            temp_res1[index] = delta
            self.breadth = 0

            temp_res1 = "".join(map(str, temp_res1))
            Ab_nt = Seq(temp_res1)
            Ab_aa = str(Ab_nt.translate())
            end1 = timer()
            mutate_time1 = end1 - start1
            start2 = timer()

            
            if len(Ab_aa) < 217:
                print("This mutation resulted in a dead B Cell due to a start/stop codon change at position: ",
                 index)
                self.residues = None
                self.E = np.inf
            else:
                self.residues[-n_clones + cln] = deepcopy(temp_res1)
                self.E[-n_clones+cln] = self.energy(Ag[0], self.cycle_number, self.residues[-n_clones+cln][0])
            
            end2 = timer()
            mutate_time2 = end2 - start2
            # ftrack.write('%d, %d,%lf,%s,%d,%s' % (-n_clones+cln,self.cycle_number, self.E[-n_clones + cln], delta, index, str(Ab_aa)))
            # ftrack.write('\n')
            # ftrack.flush()
            cln += 1


    def shm(self, Ag, cycle_number):
        print("=== Running SHM ===")
        """ Run somatic hypermutation and return self + new B cell clones. """

        # get number of cells that mutate, remove from the clone
        new_clones = []
        n_muts = np.random.binomial(self.nbs, self.p_mut)
        self.nbs -= n_muts
        
        # get number of CDR vs framework (FR) mutations
        if self.p_CDR == 1.0:
            n_CDRs = n_muts
        else:
            n_CDRs = np.random.binomial(n_muts, self.p_CDR)
        n_FR = n_muts - n_CDRs
        
        # process CDR mutations
        outcomes = self.vectorized_multinomial(ns=n_CDRs, condp=self.shm_condp)  # die, silent, affected
        n_die = outcomes[:, 0]
        n_silent = outcomes[:, 1]
        n_affect = outcomes[:, 2]
        
        self.nbs -= n_die 
        self.nbs += n_silent  # add silent mutations to the parent clone
        
        n_total = n_affect.sum() # number of new clones
        
        start = 0
        #TODO: Change this to the implementation below/amdrl for speed-up
        
        if len(self.nbs) == 1:
            new_residues = np.tile(self.residues, (n_total, 1))
        else:
            new_residues = np.tile('B' * len(self.residues[0][0]), (n_total, 1))
            start = 0            
            for res, stop in zip(self.residues, n_affect.cumsum()):
                new_residues[start:stop,:] = res
                start = stop

        # Add clones to be mutated in mutate_CDR()
        self.clone(new_residues, nbs=int(1))
        
        # Mutate clones
        self.mutate_CDR(Ag, new_residues.shape[0])
        
        # Prune clones
        self.prune()

    """
    Ideally we can use Brenden's implementation for the speed-up
    """
        # Each entry in n_affect determines how many copies of that clone (with nbs=1)
        # BP: Using np.tile was very slow. The version below achives parity and is over 10x faster.
        # n = n_affect.sum()  # Number of rows to add
        # new_residues = np.empty(shape=(n, len(self.residues)), dtype=np.float32)
        # start = 0
        # for res, stop in zip(self.residues, n_affect.cumsum()):
        #     new_residues[start:stop] = res
        #     start = stop
        #
        #     # Add the clones
        #     self.clone(new_residues, nbs=1)
        #
        #     # Mutate the newly added clones
        #     self.mutate_CDR(Ag, new_residues.shape[0])
        #
        #     # Prune clones with nbs <= 0
        #     self.prune()
    def append_mem_cells(self, new_memory_cells):
        """ 
        Function to append memory cells rather than overwrite them 

            Inputs:
                self (current memory_cells);
                new_memory_cells generated in GerminalCenter.run_recycle()
            Output:
                Appended memory_cells population to choose new founders from all memory cells
        """
        self.residues = np.vstack([self.residues, new_memory_cells.residues])
        self.nbs = np.concatenate([self.nbs, new_memory_cells.nbs])
        self.E = np.concatenate([self.E, new_memory_cells.E])
        self.ag_captured = np.concatenate([self.ag_captured, new_memory_cells.ag_captured])
        self.div = np.concatenate([self.div, new_memory_cells.div])

    def prune(self):
        print("=== Pruning ===")
        """
        Prune the dead B Cell populations
        """
        alive = self.nbs > 0
        print("B Cell Populations Pre-Pruning: ", self.nbs)
        self.ag_captured = self.ag_captured[alive]
        self.div = self.div[alive]
        self.nbs = self.nbs[alive]
        self.residues = self.residues[alive]
        self.E = self.E[alive]
        self.alive = alive
    
    # Adapted from: https://stackoverflow.com/questions/55818845/fast-vectorized-multinomial-in-python
    # NOTE: This implementation mutates ns
    def vectorized_multinomial(self, ns, condp):
        """
        Sample from the multinomial distribution with experiment sizes but the same probabilities.

        ns : experiment sizes
        condp : conditional probabilities for each outcome, pre-computed as p / p.cumsum()
        """
        if isinstance(ns, int):
            out = np.zeros(shape=(1, len(condp)), dtype=int)
        else:
            out = np.zeros(shape=(len(ns), len(condp)), dtype=int)
        for i in range(len(condp) - 1, 0, -1):
            binsample = np.random.binomial(ns, condp[i])
            out[:, i] = binsample
            ns -= binsample
        out[:, 0] = ns
        return out
    
    def selectObservations(self, observations, keys, values):
        assert len(keys) == len(values), 'Miss-matched lengths for selectObservations!'
        selected = []
        for obs in observations:
            allowed = True
            for i in range(len(keys)):
                if obs[keys[i]] != values[i]:
                    allowed = False
                    break
            if allowed: selected.append(obs)
        return selected
