# Atomistic Affinity Maturation Model — codon project

A branch of the atomistic affinity maturation model
(https://doi.org/10.1371/journal.pcbi.1009391), adapted for the codon-level work described
below. The germinal centre reaction, nucleotide-level somatic hypermutation, and S5F
targeting/substitution model all follow the published version.

## Key changes from the published model

Two substantive changes were made to the published model before this branch:

**1. Sum-of-single-points (SSP) ΔΔG scoring** (`score.py`, `singlepoint/`)

The published model scores each antibody–antigen pair structurally: build a homology model of
the complex with Modeller from a PDB template, then predict ΔG with a trained regression model
(`ppdg`). That is accurate but far too slow to run large sweeps, hence the cached model
directory and `master_list.txt`.

The `singlepoint` scoring function instead reads a precomputed table of single-point ΔΔG values
indexed by amino acid position and mutant identity, sums the ΔΔG of every mutation separating a
clone from its founder, and adds the founder's reference ΔG (`WT`). Orders of magnitude faster,
at the cost of assuming additivity across mutations (no epistasis).

**2. Explicit antigen capture** (`lymph_node.py`)

The GC now holds a finite antigen pool (`ag_start`) that depletes over the course of an
injection, rather than selecting against a fixed concentration.

- After binding selection, `assign_ag()` distributes antigen to clones in descending order of
  capture probability, capped at 200 molecules per B cell (100 BCRs × 2, after Amitai et al.,
  *PLoS Comput Biol* 2018).
- Dark-zone division rate then scales with antigen captured via `calc_expansion_rate()` (up to
  `div_max = 6`) instead of being fixed.
- Binding selection uses the *current* pool (`ag_left / ag_start`) as the effective
  concentration, and an injection ends once the pool falls below 20% of its starting size.

---

# Codon project (this branch)

## Goal

Test how the **nucleotide** sequence of the founder B cell affects affinity maturation
outcomes when the **amino acid** sequence is held fixed — that is, synonymous codon variants
of the same germline antibody.

The model separates the two levels cleanly, which is what makes the question tractable:

- Somatic hypermutation acts on nucleotides, and mutation sites are drawn from the S5F 5-mer
  mutability/substitution tables (`Mutability.csv`, `Substitution.csv`). Codon choice changes
  *which mutations are reachable*.
- With `"scoring_function": "singlepoint"`, binding energy depends only on the amino acid
  sequence. Synonymous variants therefore share an identical fitness landscape and identical
  starting ΔG.

So codon choice influences the outcome only through mutational accessibility, with no
confounding change to the energy function. (Using `ppdg` instead would add structural-model
noise on top and blur this separation.)

## Design notes for building codon variants

- **Only CDR nucleotides ever mutate.** `BCell.mutate_CDR` collects 5-mers whose *center* falls
  inside a CDR, and `n_FR` is computed in `shm` but never applied to the sequence. Framework
  synonymous changes are inert — except within 2 nt of a CDR boundary, where they still shift
  the 5-mer context of CDR-edge positions.
- **The last 2 nt of the sequence can never mutate**, since the 5-mer scan stops at
  `i < len(seq) - 4`. This is the tail of `CDRL3_nt`.
- **Region lengths must be preserved** across variants. Lengths are read from the first
  antibody in the config and used to locate every CDR boundary.
- **Founder names must end in `_<integer>`** (e.g. `founder_1`). `lymph_node.py` parses the
  trailing integer as the `parent_id` used for per-founder ΔΔG lookup.
- For head-to-head runs, set `multiple_founders: true` and pass the *same* ΔΔG file once per
  founder in `singlepoint_file`; `parent_id` is tracked through cloning, mutation, and the
  HDF5 log, so lineages can be attributed to a variant.

## Seeding and reproducibility

Runs are stochastic by default and are meant to be — replicates are the point. Seeds are
recorded rather than fixed, so that any *individual* run can be revisited later without giving
up run-to-run variation.

By default a fresh seed is drawn from OS entropy at the start of each run and is:

- printed to stdout (so it lands in the SLURM `.out` file), and
- written into `config['experiment']['seed']`, which is saved to `config.json` in that run's
  log directory.

To repeat a specific run, supply its recorded seed:

```bash
python -u -m run --config=config/my_config.json --seed 3857938007
```

or add `"seed": 3857938007` to the `experiment` block. Precedence is
`--seed` > config `"seed"` > freshly drawn. A single seed covers the whole run, all GCs
included, since they share the global RNG.

Note that this reproduces runs using `singlepoint` scoring. The `ppdg` path shells out to
Modeller and reads a shared model cache, so it is not reproducible from the seed alone.

## Changes from the starting branch

### 1. `mutate_CDR`: rebuild the 5-mer pool for each clone (`lymph_node.py`)

`Fivemer` and `Fivemer_index` were initialized once *before* the per-clone loop and never
cleared inside it, so candidate 5-mers accumulated across clones. When mutating clone *k*, the
sampling pool held all previous clones' 5-mers alongside its own, so the chosen mutation site
and substitution context could come from a different clone's sequence. Since all sequences are
the same length the position index stayed valid, producing incorrect mutations rather than an
error.

The effect is silent when every clone in a batch shares a sequence (repeated 5-mers normalize
out) and grows with clonal diversity. Under `multiple_founders`, batches are heterogeneous from
the first cycle.

Both lists are now initialized inside the loop.

> Note: this behavior was inherited from `mixture.py`, not introduced by the refactor. Results
> generated before this fix are affected.

### 2. `shm`: stop subtracting lethal CDR mutants twice (`lymph_node.py`)

`self.nbs -= n_muts` already removes every mutating cell from its parent clone. The following
`self.nbs -= n_die` removed the lethal fraction a second time, from cells that had not mutated.
For a clone of 100 with 10 mutating and 3 lethal, the population became 94 instead of 97.

That line is removed. The accounting is now: mutators leave the clone, silent mutants return,
affinity-affecting mutants become new clones, and lethal (`n_die`) plus framework (`n_FR`)
mutants are the deaths.

> Note: this line is **not** present in `mixture.py` and was introduced during the refactor.
> Its intent is unconfirmed. If earlier results from this lab were generated with it, they are
> not directly comparable to runs made after this change.

### 3. Seed plumbing (`run.py`, `lymph_node.py`, `config/generic_config.json`)

Added `resolve_seed()` and a `--seed` option to `run.py`, and a `"seed"` field to the generic
config (default `null`). The bare `np.random.seed()` at import time in `lymph_node.py` is
removed, since seeding now happens once in `run.py` where the value can be recorded. Both
`numpy` and the stdlib `random` module are seeded — `random` is used for epitope selection in
`BCell.energy`. Behavior is unchanged when no seed is given, other than being recorded.

## Known issues, not addressed

- `Mutability.csv` and `Substitution.csv` are re-read from disk and scanned linearly for every
  5-mer of every clone. Correct but slow; loading them once into a dict keyed by 5-mer would
  speed up large sweeps.
- `pick_memCells` uses `np.int`, removed in NumPy >= 1.24. Reached only on multi-injection runs.
- `config/__init__.py` uses `collections.Mapping`, removed in Python 3.10.
