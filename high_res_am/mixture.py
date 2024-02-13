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

warnings.filterwarnings("ignore")

np.random.seed()  # Set NumPy random seed

# Set to True to turn off behavior changes since BCell refactor
# Set to False to turn on new behavior changes:
#   1. Memory cells also added on last cycle.
#   2. Memory founders proliferate to 512 instantly.
#   3. Memory founders are selected at the clonal level (not weighted by nbs), without replacement.
REVERT = False

###### Global parameters ######

p_mut = 0.14  # probability of mutation per division round
p_CDR = 0.85  # probability of mutation in the CDR region
p_CDR_lethal = 0.30  # probability that a CDR mutation is lethal
p_CDR_silent = 0.00  # probability that a CDR mutation is silent
p_CDR_affect = 1. - p_CDR_lethal - p_CDR_silent  # probability that a CDR mutation affects affinity
p_recycle = 0.70  # probability that a B cell is recycled
p_exit = 1. - p_recycle  # probability that a B cell exits the GC
nb_GC_founders = [1, 2, 3]
nb_GC_founders_I1 = 1
nb_GC_founders_I2or3 = 3
epsilon = 1e-16
activation_energy = -9.011959719155849
nb_trial = 1
energy_scale = 3.5 / 0.6  # inverse temperature : 1/(8.3145/1000*310)
help_cutoff = 0.70  # only B cells in the top (help_cutoff) fraction of binders receive T cell help

# Upload dictionary of antigens
import dictionary_little_code

importlib.reload(dictionary_little_code)
from dictionary_little_code import dicAgs, dicconc, dicGCDur, flag, time_Steps


print('flag = ', flag)

# Upload seeding cells
import seedingBcell

importlib.reload(seedingBcell)
from seedingBcell import seedingCells, VRC01GL_nt, VRC01GL_nt_str

# Upload founder cells
if flag == 2 or flag == 3:
    import founders

    importlib.reload(founders)
    from founders import founderseq

lenFRWH1 = len(VRC01GL_nt[0][0])
lenCDRH1 = len(VRC01GL_nt[0][1])
lenFRWH2 = len(VRC01GL_nt[0][2])
lenCDRH2 = len(VRC01GL_nt[0][3])
lenFRWH3 = len(VRC01GL_nt[0][4])
lenCDRH3 = len(VRC01GL_nt[0][5])
lenFRWL1 = len(VRC01GL_nt[1][0])
lenCDRL1 = len(VRC01GL_nt[1][1])
lenFRWL2 = len(VRC01GL_nt[1][2])
lenCDRL2 = len(VRC01GL_nt[1][3])
lenFRWL3 = len(VRC01GL_nt[1][4])
lenCDRL3 = len(VRC01GL_nt[1][5])

# ftrack = open('output-track.txt', 'w')
# ftrack2 = open('output-track2.txt', 'w')
# ftrack3 = open('output-track3.txt', 'w')
# ftrack4 = open('output-track4.txt', 'w')
# ftotTEST = open('output-totalTEST.csv', 'w')
# fseed = open('seed.csv', 'w')


# Adapted from: https://stackoverflow.com/questions/55818845/fast-vectorized-multinomial-in-python
# NOTE: This implementation mutates ns
def vectorized_multinomial(ns, condp):
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


# Conditional p values used in vectorized_multinomial() are repeated often, so compute them once here.
p = np.array([p_CDR_lethal, p_CDR_silent, p_CDR_affect])
ps = p.cumsum(axis=-1)
with np.errstate(divide='ignore', invalid='ignore'):
    condp = p / ps
condp[np.isnan(condp)] = 0.0
shm_condp = condp


class BCell:
    # !! FIX !! #
    def __init__(self, residues, nbs, **kwargs):
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
        # NOTE: this is the number of mutations in the GL CODON sequence,
        # not AA sequence (# aa is < or = # nt mutations)
        self.residues = residues
        self.nbs = nbs
        if 'nb_CDR_mut' in kwargs:
            self.nb_CDR_mut = kwargs['nb_CDR_mut']
        else:
            self.nb_CDR_mut = 0

        if 'residues' in kwargs:
            self.residues = residues
        # else:
        #     print('residues not recognized as an input argument')

        if 'cycle_number' in kwargs:
            self.cycle_number = kwargs['cycle_number']
        else:
            self.cycle_number = 2

        if 'antigens' in kwargs:
            self.antigens = np.array(kwargs['antigens'])
        else:
            self.antigens = np.array(dicAgs[self.cycle_number])

        if 'E' in kwargs:
            self.E = kwargs['E']
        else:
            if isinstance(nbs, int) and isinstance(residues, str):
                self.residues = str(self.residues)
                self.E = self.energy(self.antigens[0], self.cycle_number, self.residues)
            else:
                self.E = np.inf
            # assuming single Ag case!!!!### ###FIX###

        if 'breadth' in kwargs:
            self.breadth = np.array(kwargs['breadth'])
        else:
            self.breadth = 0

        if 'nb_Ag' in kwargs:
            self.nb_Ag = kwargs['nb_Ag']
        else:
            self.nb_Ag = np.shape(self.antigens)[0]

        # if 'last_bound' in kwargs:
        #     self.last_bound = kwargs['last_bound']
        # else:
        #     self.last_bound = np.random.multinomial(self.nbs, pvals=[1 / float(self.nb_Ag)] * self.nb_Ag)

        if 'mut_res_id' in kwargs:
            self.mut_res_id = kwargs['mut_res_id']
        else:
            self.mut_res_id = 0

        if 'delta_res' in kwargs:
            self.delta_res = kwargs['delta_res']
        else:
            self.delta_res = 0

        if 'generation' in kwargs:
            self.generation = kwargs['generation']
        else:
            self.generation = 0

        if 'history' in kwargs:
            self.history = kwargs['history']
        else:
            self.history = {'generation': [self.generation], 'cycle_number': [self.cycle_number],
                            'residues': [self.residues], 'nb_CDR_mut': [self.nb_CDR_mut],
                            'mut_res_id': [self.mut_res_id], 'E': [self.E], 'delta_res': [self.delta_res]}

        if residues is None:
            residues = seedingCells[0]
        self.residues = residues
        self.nbs = nbs
        self.E = np.array(self.E).reshape((1,))


    def clone(self, residues, nbs):
        """ Add clones with given residues and nbs """
        print('starting clone()')
        # If nbs is an int, repeat that value for all new clones
        if isinstance(nbs, int):
            nbs = np.repeat(nbs, residues.shape[0])
            E = np.repeat(self.E[0], residues.shape[0])

        self.residues = np.vstack([self.residues, residues])
        self.nbs = np.concatenate([self.nbs, nbs])
        self.E = np.concatenate([self.E, E])
        return self

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
        print('starting energy()')
        # Convert to Bio nt_seq
        Ab_nt = Seq(residues)

        # Convert to AA
        Ab_aa = str(Ab_nt.translate())
        
        # Format for chains in Modeller
        Ab = Ab_aa[0:118] + '/' + Ab_aa[118:]
        print(self.nbs, cycle_number, self.nb_CDR_mut)
        if "*" not in Ab_aa:
            print('Ag:', Ag)
            # print('Ab=', Ab)
            print("Ab:", Ab)
            dg = compute_descriptors_kayla.main(Ag, Ab)
        else:
            dg = 100
        return dg

    
    def divide(self, cycle_number):
        """ Run one round of division. """
        print('starting divide()')
        self.nbs *= 2
        self.generation += 1
        self.cycle_number += 1


    def calculate_breadth(self, testpanel, threshold, panelSize):
        test_energies = np.inner(testpanel, self.residues)
        return np.mean(test_energies > threshold, axis=0)

    def mutate_CDR(self, Ag, n_clones):

        """
        Change in energy due to affinity-affecting CDR mutation. Only one residue mutates.

        Targeting model: choose nt to mutate based on 5mer mutability and sequence occurrences

        Substitution model: choose nt to mutate to based on probability of middle nt in 5mer
                            mutating to any other nt
        """
        print("Starting mutate_CDR()")
        start1 = timer()
        global picked_list2
        Fivemer = []
        Fivemer_index = []
        
        residues = self.residues[-n_clones:]
        print('n_clones:', n_clones)
        cln = 0
        for clone_residues in residues:
            clone_residues = clone_residues[0]
            for i in range(len(str(clone_residues))):  # Iterating over the entire sequence
                if i < (
                        len(str(clone_residues)) - 4):
                    # If i is less than four spaces from the end of the sequence
                    # (i.e., there is a fivemer available)
                    if lenFRWH1 <= (i + 2) < lenFRWH1 + lenCDRH1:
                        # If the length of the first heavy FRW region is less than or equal to i + 2,
                        # and if i+2 is less than the combined length of FRWH1 and CDRH1.
                        # So essentially, if there is still a Fivemer in the FRWH1,
                        # and i is within the bounds of the first mutation region
                        Fivemer_index.append(i + 2)
                        Fivemer.append(
                            str(clone_residues)[i + 0:i + 5])
                        # Appending the section of the sequence to check for mutability
                    elif lenFRWH1 + lenCDRH1 + lenFRWH2 <= (i + 2) < lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 <= (
                            i + 2) < lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 <= (
                            i + 2) < lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 + lenCDRL1:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 + lenCDRL1 + lenFRWL2 <= (
                            i + 2) < lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 + lenCDRL1 + lenFRWL2 + lenCDRL2:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])
                    elif lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 + lenCDRL1 + lenFRWL2 + lenCDRL2 + lenFRWL3 <= (
                            i + 2) < lenFRWH1 + lenCDRH1 + lenFRWH2 + lenCDRH2 + lenFRWH3 + lenCDRH3 + lenFRWL1 + lenCDRL1 + lenFRWL2 + lenCDRL2 + lenFRWL3 + lenCDRL3:
                        Fivemer_index.append(i + 2)
                        Fivemer.append(str(clone_residues)[i + 0:i + 5])

            reader = csv.DictReader(open('Mutability.csv'))
            mutTable = [obs for obs in reader]
            a = [float(o['Mutability']) for o in mutTable]
            mutSum = np.sum(a)

            probs = []
            for i in Fivemer:
                FivemerMutTable = selectObservations(mutTable, ['Fivemer'], [str(i)])
                a = [float(o['Mutability']) for o in FivemerMutTable]
                probs.append(a[0])
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
                FivemerMutTable2 = selectObservations(mutTable2, ['Fivemer'], [str(i)])
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
            print("index/delta:", index, delta)
            
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
            ftrack.write('%d, %d,%lf,%s,%d,%s' % (-n_clones+cln,self.cycle_number, self.E[-n_clones + cln], delta, index, str(Ab_aa)))
            ftrack.write('\n')
            ftrack.flush()
            cln += 1

    def shm_old(self, Ag, cycle_number):
        """ Run somatic hypermutation and return self + new B cell clones. """
        print("Starting shm()")
        # get number of cells that mutate, remove from the clone
        new_clones = []

        n_muts = np.random.binomial(self.nbs, p_mut)
        self.nbs -= n_muts

        # get number of CDR vs framework (FR) mutations
        if p_CDR == 1.0:
            n_CDRs = n_muts
        else:
            n_CDRs = np.random.binomial(n_muts, p_CDR)
        n_FR = n_muts - n_CDRs

        # process CDR mutations
        outcomes = vectorized_multinomial(ns=n_CDRs, condp=shm_condp)  # die, silent, affected
        n_silent = outcomes[:, 1][0]
        n_affect = outcomes[:, 2][0]

        self.nbs += n_silent  # add silent mutations to the parent clone

        print("naffect= ", n_affect)
        # ftrack2.write('n_affect = {:d},{:d},{:d}'.format(n_affect, cycle_number, self.nbs))
        # ftrack2.write('\n')
        # ftrack2.flush()

        for j in range(len(self.nbs)):
            ftrack.write('clone: %d' % j)
            ftrack.write('\n')
            ftrack.flush()
            print("Clone: ", j)
            for i in range(n_affect):
                b = self.clone(self.residues[j], self.nbs[j])
                Ag = b.antigens[0]

                b.mutate_CDR(Ag, cycle_number)
                print("mutate_cdr b.nbs:", b.nbs)
                if b.residues is not None and b.nbs > 0:
                    new_clones.append(b)
                    self.residues[j] = b.residues
                    self.nbs[j] = b.nbs
                    self.nb_CDR_mut[j] += 1
                    self.E[j] = b.E
                    print("energy in SHM(): ", self.E[j])

        self.prune()

        return self

    def shm(self, Ag, cycle_number):
        """ Run somatic hypermutation and return self + new B cell clones. """
        print("Starting shm()")
        # get number of cells that mutate, remove from the clone
        new_clones = []
        print(self.residues)
        # Should we use the vectorized multinomial?
        n_muts = np.random.binomial(self.nbs, p_mut)
        self.nbs -= n_muts

        # get number of CDR vs framework (FR) mutations
        if p_CDR == 1.0:
            n_CDRs = n_muts
        else:
            n_CDRs = np.random.binomial(n_muts, p_CDR)
        n_FR = n_muts - n_CDRs

        # process CDR mutations
        outcomes = vectorized_multinomial(ns=n_CDRs, condp=shm_condp)  # die, silent, affected
        n_silent = outcomes[:, 1]
        n_affect = outcomes[:, 2]

        self.nbs += n_silent  # add silent mutations to the parent clone

        n_total = n_affect.sum() # number of new clones
        print('n_total =', n_total)
        print("n_affect =", n_affect)
        start = 0

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
        for item in self.nbs:
            ftrack.write("nbs: %d" % item)
        ftrack.write('\n')
        ftrack.flush()

        # Prune clones that are dead
        self.prune()

        return self

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

    def prune(self):
        """
        Prune the dead B Cell populations
        """
        alive = self.nbs > 0
        self.nbs = self.nbs[alive]
        self.residues = self.residues[alive]
        self.E = self.E[alive]


class GerminalCenter(object):

    def __init__(self, body, debug=0):
        self.debug = debug
        self.memory_cells = BCell(residues=None, nbs=1)
        self.memory_founders = BCell(residues=None, nbs=1)
        # Single Antigen
        # if flag == 1:

        # Get the seeding cells from a random choice w/o replacement
        id_seeding_cells = np.random.choice(len(seedingCells), nb_GC_founders_I1, replace=False)

        # Residues of seeding B Cells from germ lines
        # residues = seedingCells[id_seeding_cells[0]]
        residues = seedingCells[0]

            # self.memory_founders = BCell(residues=None)

        # Multiple Antigens (2 or 3)
        # elif flag == 2 or flag == 3:

            # Get founders for each GC from a random choice
            # Aren't all of these the same??
            # id_founding_cells = np.random.choice(len(founderseq), nb_GC_founders_I2or3, replace=False)

            # Initialize mutation counter
            # cdr_mut = []
            # Iterate over each of the sequences in the founderseq set
            # for seq in founderseq:

                # Tally the changes between the current nt_seq and the GL_nt_seq
                # mutations_aa = [i for i in range(len(VRC01GL_nt_str)) if VRC01GL_nt_str[i] != seq[i]]
                # cdr_mut.append(len(mutations_aa))

                # Generate memory_founders with the current sequences
                # Is this just to keep track of the changes?
            # self.memory_founders = [BCell(residues=founderseq[id_founding_cells[i]], nb_CDR_mut=cdr_mut[id_founding_cells[i]]) for i in
            #            range(nb_GC_founders_I2or3)]

        # print("The residues in GC __init__ are:", residues)
        self.B_cells = BCell(residues=residues, nbs=128)
        # self.memory_cells = BCell(residues=None)
        # self.plasma_cells = BCell(residues=None)

        # all founders replicate and are recycled
        # maximum number of cells in the GC (= initial population size)
        # self.GC_size_max = self.B_cells.nbs
        if isinstance(self.B_cells.nbs, int):
            self.GC_size_max = self.B_cells.nbs
        else:
            self.GC_size_max = self.B_cells.nbs.sum()
        self.cycle_number = 2
        self.nb_cycle_max = 99
        self.n_injections = 0  # Number of injections so far

    def injection(self, dicAgs, dicAbs):
        self.n_injections += 1

        done = False

        # Set the antigens and concentrations
        self.Ag = dicAgs[0]['seq']
        self.conc = dicAgs[0]['conc']

        # Second injection onward, use memory founders
        if self.n_injections > 1:
            self.cycle_number = 0  # Reset cycle number
            self.plasma_cells = BCell(residues=None)  # Reset plasma cells

            # Set current B Cells to be memory founders
            self.B_cells = self.memory_founders

            # If there aren't any B Cells, exit since the GC has expired
            if len(self.B_cells.nbs) == 0:
                if self.debug > 0:
                    print("GC failed; can't perform injection. Returning True")
                return True

        # The value in this while loop should be a global var that is set at runtime
        while self.cycle_number < self.nb_cycle_max:
            self.run_GC_cycle()  # Run the GC cycle once
            print(self.cycle_number, self.B_cells.nbs.sum())

            """
            test_state() is checking to see if the current GC has exceeded the maximum
            size. It will also check if the maxmimum cycle_number has been reached. If
            neither flag is raised, it will then add 1 to the current cycle_number.
            """
            done, exceed_size_flag = self.test_state()
            if done:
                if self.debug > 1:
                    print("Finished after cycle {}".format(self.cycle_number))
            elif exceed_size_flag:
                self.memory_founders = self.pick_memCells()
            break

        return done

    def run_GC_cycle(self):
        """ Run one cycle of the GC reaction. """
        print("starting run_GC_cycle()")
        print('nbs: ', self.B_cells.nbs)
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

    def run_dark_zone(self, nb_rounds=2):
        """ B cells proliferate and undergo SHM in the dark zone. """
        print('starting run_dark_zone()')
        for i in range(nb_rounds):
            self.B_cells.divide(self.cycle_number)
            self.B_cells.shm(self.Ag, self.cycle_number)

    def run_binding_selection(self):
        """ Select B cells for binding to antigen. """
        print('Starting run_binding_selection()')
        for i in range(len(self.B_cells.nbs)):
            Ag_bounds = np.exp(
                energy_scale * (self.B_cells.E[i] - activation_energy))  # Shape: (n_clones,)
            factors = self.conc * Ag_bounds
            langmuir_conjs = 1. / (1. + factors)
            n_dies = np.random.binomial(self.B_cells.nbs[i], langmuir_conjs)
            print("B-Cells dying in binding selection:", n_dies)
            self.B_cells.nbs[i] -= n_dies
        self.B_cells.prune()

    def run_help_selection(self):
        """ Select B cells to receive T cell help. """
        print('Starting run_help_selection()')
        binding_energies = self.B_cells.E
        
        # for i in range(len(self.B_cells.nbs)):
        #     binding_energies.append(self.B_cells.energy(self.Ag[0], self.cycle_number, self.B_cells.residues[i]))
        binding_energies = np.array(binding_energies)
        sort = np.argsort(binding_energies)[::-1]

        sorted_binding_energies = binding_energies[sort]
        sorted_nbs = self.B_cells.nbs[sort]
        nb_tot = self.B_cells.nbs.sum()
        if nb_tot > 0:
            cut_idx = int(help_cutoff * nb_tot) - 1
            energy_idx = np.argmax(np.cumsum(sorted_nbs) > cut_idx)
            # Returns index of first instance in which cumulative sum of nbs is greater than the cutoff value
            energy_cutoff = sorted_binding_energies[energy_idx]

            # Kill all B cells below threshold
            dies = binding_energies < energy_cutoff
            self.B_cells.nbs[dies] = 0

            # Kill a fraction of B cells exactly at threshold
            # Tie deaths are split evenly across all clones exactly at threshold
            ties = binding_energies == energy_cutoff
            n_die_tie = self.B_cells.nbs.sum() - cut_idx
            # Total number of cells to remove for clones exactly at threshold
            fraction_die_tie = n_die_tie / self.B_cells.nbs[ties].sum()
            # Fraction of cells to remove for each clone exactly at threshold
            self.B_cells.nbs[ties] = np.multiply(self.B_cells.nbs[ties], (1 - fraction_die_tie), casting="unsafe")
            # Note: Can't use *= because we're multiplying an int by a float
            self.B_cells.prune()

    def run_recycle(self):
        """ Randomly select B cells to be recycled back into the GC or to become memory cells. """
        print("starting run_recycle()")
        n_exits = np.random.binomial(self.B_cells.nbs, p_exit)
        self.B_cells.nbs -= n_exits
        new_memory_cells = BCell(residues=self.B_cells.residues, nbs=n_exits)
        new_memory_cells.E = self.B_cells.E
        self.B_cells.prune()

        if REVERT:  # Old behavior: Don't add memory cells for the last cycle
            GC_size = self.B_cells.nbs.sum()
            if (self.cycle_number == self.nb_cycle_max - 1) or (GC_size > self.GC_size_max):
                pass
            else:
                self.memory_cells.clone(new_memory_cells.residues, new_memory_cells.nbs)  # Always add memory cells
        else:
            self.memory_cells = new_memory_cells  # Always add memory cells

    def pick_memCells(self):

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

        if flag == 1:
            i = 0
        else:
            i = 2

        if REVERT:
            if n_mem_clones < nb_GC_founders[i]:
                if self.debug > 0:
                    print("Not enough memory cells to form founders!")
                return BCell(residues=None, nbs=0)
            p = None
            replace = False
        else:
            if nb_tot < nb_GC_founders[i]:
                if self.debug > 0:
                    print("Not enough memory cells to form founders!")
                return BCell(residues=None, nbs=0)
            p = nbs / nb_tot # Probabilities for choosing each clone (proportaional to nbs)
            replace = True # Approximation

        id_new_founders = np.random.choice(n_mem_clones, nb_GC_founders[i], p=p, replace=replace)
        residues = self.memory_cells.residues[id_new_founders]

        if REVERT:
            nbs = self.memory_cells.nbs[id_new_founders]
        else:
            nbs = 128  # Instantly proliferate to 128

        new_founders = BCell(residues=residues, nbs=nbs)
        return new_founders

    def test_state(self):
        GC_size = self.B_cells.nbs.sum()  # total number of cells in the GC
        if (self.cycle_number == self.nb_cycle_max - 1) or (
                GC_size > self.GC_size_max):  # at the end, all B cells exit the GC
            self.plasma_cells = self.B_cells

        if (GC_size <= 0):
            exceed_size_flag = False
            done = True
            if self.debug > 0:
                print("Done: GC size = 0 (all dead)")
        elif (GC_size > self.GC_size_max):
            exceed_size_flag = True
            done = False
        else:
            exceed_size_flag = False
            done = False

        self.cycle_number += 1

        return done, exceed_size_flag


class AffinityMaturation(object):

    def __init__(self, n_GCs, debug=0):
        self.gcs = [GerminalCenter(body=self, debug=debug) for _ in range(n_GCs)]

    def injection(self, Ag, conc):
        # Call injection for each GC
        dones = [gc.injection(Ag, conc) for gc in self.gcs]

        return all(dones)


###### Main functions ######

def main(verbose=False):

    # Run multiple trials and save all data to file
    start = timer()
    n_GCs = 1
    n_injections = 1
    Ag = dicAgs
    conc = dicconc
    kwargs = {
        ''
    }
    simulation = AffinityMaturation(n_GCs=n_GCs, debug=0)
    for i in range(n_injections):
        print("Starting injection {} of {}".format(i+1, n_injections))
        done = simulation.injection(Ag[time_Steps[i]], conc[time_Steps[i]])
        nbs = np.sum([gc.B_cells.nbs.sum() for gc in simulation.gcs])
        print('Total B Cells:', nbs)
        if done:
            break
    end = timer()
    print('\nTotal time: %lfs, average per GC %lfs' % ((end - start), (end - start) / float(n_GCs)))


def selectObservations(observations, keys, values):
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


def find_nb_Ag(antigens):
    nb_Ag = np.shape(antigens)[0]
    return nb_Ag


def pick_memCells_for_new_GC(memory_cells, nb_GC_founders):
    n_mem_cells = len(memory_cells)
    id_new_founders = np.random.choice(n_mem_cells, nb_GC_founders, replace=False)
    new_founders = [memory_cells[id_new_founders[i]] for i in range(nb_GC_founders)]
    return new_founders


def run_breadth_calculation(panel_energies, threshold, panelSize):
    average = np.mean(panel_energies)
    variance = np.var(panel_energies)
    breadth = float(sum(x > threshold for x in panel_energies)) / panelSize
    return average, variance, breadth


if __name__ == '__main__': main()
