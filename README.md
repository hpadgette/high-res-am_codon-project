# Atomistic Affinity Maturation Model — codon project

A branch of the atomistic affinity maturation model
(https://doi.org/10.1371/journal.pcbi.1009391), adapted for the codon-level work described
below. The germinal centre reaction, nucleotide-level somatic hypermutation, and S5F
targeting/substitution model all follow the published version.

## Key changes from the published model

Two substantive changes are made to the published model with this branch:

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

## Where things live

| To change... | Edit |
|---|---|
| Antigen, founders, any model parameter | your experiment JSON in `config/` |
| Parameter defaults | `config/generic_config.json`, `config/generic_hiv_no_comments.json` |
| The GC reaction (selection, mutation, recycling) | `lymph_node.py` |
| How binding energy is computed | `score.py` |
| Experiment setup / how founders are assembled | `affinity.py` |
| CLI, seeding, log directory layout | `run.py` |

Each module has a docstring at the top explaining its role and what a new user would change
there. Run from the `high_res_am/` package directory — `Mutability.csv` and
`Substitution.csv` are opened by relative path.

`json.load` cannot parse `//` comments, so **any config that is actually loaded must be plain
JSON**. `generic_hiv.json` is the commented reference; `generic_hiv_no_comments.json` is the
one the loader reads. Edit one, edit both.

## Configuration reference

### `experiment`

| Key | Meaning |
|---|---|
| `logging` | Write log files at all |
| `log_level` | Which logs to write: `"injection"`, `"cycles"`, `"B Cell"` |
| `exp_name` | Name of the run's log directory (a timestamp is appended) |
| `num_GCs` | Number of independent germinal centres |
| `num_injections` | Prime/boost count — **see warning below** |
| `debug` | `0` for normal runs; higher values skip log setup and print more |
| `seed` | RNG seed. `null` draws a fresh one and records it |
| `score_config.scoring_function` | `"singlepoint"`, `"ppdg"`, or `"ppdx"` (unimplemented) |
| `score_config.singlepoint_file` | Path to a ΔΔG CSV, or a list of them for multiple founders |
| `score_config.WT` | Founder reference ΔG, or a list paired positionally with the files |

### `antigen` / `antibody`

| Key | Meaning |
|---|---|
| `antigen_names`, `antigen_seqs` | Antigen label and amino acid sequence |
| `antibody_names` | Founder labels. Must end in `_<integer>` — that integer is the `parent_id` used for per-founder ΔΔG lookup |
| `antibody_seq` | Per-founder dict of the 12 nucleotide regions (`FRWH1_nt` … `CDRL3_nt`) |

### `maturation_parameters`

| Key | Typical | Meaning |
|---|---|---|
| `p_mut` | 0.14 | Probability a cell mutates, per dark zone round |
| `p_CDR` | 0.85 | Fraction of mutations landing in a CDR; the rest are framework hits and are treated as deaths |
| `p_CDR_lethal` | 0.30 | Of CDR hits, the fraction that kill the cell |
| `p_CDR_silent` | 0.00 | Of CDR hits, the fraction with no effect |
| `p_recycle` | 0.70 | Probability a survivor re-enters the GC rather than exiting as a memory cell |
| `p_scale` | 0.25 | Scales how much antigen a good binder captures per cycle |
| `energy_scale` | 0.9 | Selection sharpness — how strongly ΔG differences drive survival |
| `help_cutoff` | 0.70 | Fraction of cells surviving T cell help. Higher = weaker selection |
| `pop_size` | 100–200 | Starting B cell population per GC |
| `ag_start` | 50k–100k | Size of the antigen pool at the start of an injection |
| `num_GC_founders` | 10 | Founders seeding each GC |
| `multiple_founders` | **true** | **Required** — seeds the GC founder population *and* computes each founder's starting ΔG from its sequence. Set `true` even for one founder; see warning below |
| `split_founders` | false | **Required** — split the starting population across `num_GC_founders` copies of a *single* sequence. `false` unless you want that |
| `prob_metric` | 0 or list | **Required** — the run fails with a `KeyError` if it's missing. Epitope weights for multi-epitope binding, e.g. `[0.6, 0.3, 0.1]`. Use `0` for a single Ag–Ab epitope pair (any non-list means epitope 0). Only consulted by `ppdg`; inert under `singlepoint` |
| `revert` | false | Legacy memory-cell sampling. Keep `false` |

`p_CDR_affect` (= 1 − `p_CDR_silent` − `p_CDR_lethal`) and `p_exit` (= 1 − `p_recycle`) are
derived in `config/__init__.py`. Don't set them by hand.

> **Set `multiple_founders: true` even for a single founder** (with `num_GC_founders: 1`).
> The name is misleading: the flag gates the block in `BCell.__init__` that seeds the founder
> population *and* computes each founder's ΔG from its sequence. With `false`, that block is
> skipped and `__init__` falls through to `self.E = np.inf` — the founder enters the GC unable
> to bind antigen, captures nothing, and dies within a cycle or two.
>
> The first symptom looks unrelated: `ValueError: zero-dimensional arrays cannot be
> concatenated` in `append_mem_cells`, because the `memory_cells` placeholder's scalar `nbs`
> is never converted to an array. That crash is a *symptom* — fix the flag, not `lymph_node.py`.
>
> Neither `multiple_founders` nor `split_founders` has a default in `generic_config.json`, so
> both must appear in your experiment config or the run fails with a `KeyError`.

### Keys that look live but are not

These appear in the example configs and have **no effect** on a run:

| Key | Why |
|---|---|
| `activation_energy` | Overridden in `GerminalCenter.__init__` as (weakest founder ΔG) + 0.75 kcal/mol. It tracks `WT` automatically, so leave the key out rather than setting a stale number. To change selection stringency, edit that line — a larger offset is more permissive |
| `concentration` | Only read by the inactive `run_binding_selection_old()`. The live path uses `ag_left / ag_start` |
| `flag` | Sets a local variable in `pick_memCells()` that is never read |
| `cocktail` | Never read from the config |
| `num_antigens`, `delayed_release`, `time_steps`, `epsilon`, `nb_trial` | Not referenced in any active module |

> **`num_injections > 1` does not currently work.** Both calls to `pick_memCells()` in
> `GerminalCenter.injection()` are commented out (`"TODO: Change this back; only temp"`), so
> `memory_founders` is never populated and injection 2 starts with `residues=None`. Restore
> those calls before attempting prime-boost runs. Single-injection runs are unaffected.

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

## The CLK19gl → eOD-GT8 run (`codon-project/`)

`codon-project/` holds this project's experiment config and the data it was derived from.
Submit it from the inner `high_res_am/` directory, the same place `run.py` lives:

```bash
python -u -m run --config=codon-project/config_CLK19gl.json [--seed N]
```

| File | What it is |
|---|---|
| `config_CLK19gl.json` | The experiment config |
| `CLK19_final_nt_sequences.fasta` | Source nt sequences (`gl`, plus `min`/`max` synonymous recodes) |
| `CLK19_antibody_seq_derived.json` | The same sequences already split into the 12 config regions |
| `tpl_CLK19gl_eOD-GT8_modeller.pdb` | Structure the ΔΔG table was computed on |
| `CLK19_pipeline.ipynb`, `CLK19_segmentation_derivation.ipynb` | How the regions and segmentation were derived |

Key settings and why:

| Setting | Value | Note |
|---|---|---|
| `scoring_function` | `singlepoint` | AA-only energy, so synonymous variants share one landscape |
| `singlepoint_file` | `singlepoint/CLK19gl_eOD-GT8_ddg.csv` | Relative, so it resolves identically locally and on the HPC |
| `WT` | −12.257285857077804 | Germline reference ΔG, from `descriptors-med.json`. The CSV is ΔΔG *relative to this* |
| `antibody_names` | `["CLK19gl_0"]` | The `_0` suffix is mandatory — it's parsed as `parent_id` |
| `num_GC_founders` / `multiple_founders` | 1 / `true` | One founder, seeded properly. See the warning in the config reference |
| `num_injections` | 1 | Prime-boost does not currently work |

Derived values worth knowing: E_a = `WT` + 0.75 = **−11.507**, and the founder enters at
−12.257. A healthy first cycle keeps the germline lineage alive, so `worst_E` in
`cycle_logs/` should still equal `WT` early on; if the founder vanishes and `n_clones`
collapses in cycle 1, check `multiple_founders` first.

### Sequence ↔ structure mapping

The AM reading frame is **HC 360 nt / 120 aa** and **LC 300 nt / 100 aa**, mapping 1:1 onto
`tpl_CLK19gl_eOD-GT8_modeller.pdb` chain B 178–297 and chain C 401–500. The ΔΔG CSV covers
224 positions; lookups use 0–219 over the HC+LC concatenation (`__call__` strips the `/`
chain separator first), so positions 220–223 — the Cκ elbow past the reading frame — are
never read. Verified by checking that the residue absent from each position's 19 mutations is
exactly the WT residue at that position; a frame shift would break that pattern, which makes
it a quick way to re-validate any new ΔΔG table against its sequence.

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

### 4. `run_help_selection`: invert the cut so the *best* `help_cutoff` fraction survives (`lymph_node.py`)

`help_cutoff` is documented as the fraction of cells surviving T cell help, but the cut was
computed from the wrong end of the distribution. `np.argsort(binding_energies)[::-1]` orders
clones by *descending* ΔG — worst binder first — so walking the cumulative sum to
`help_cutoff * nb_tot` lands the threshold at the 70th percentile counting up from the bottom.
The result was that `help_cutoff = 0.70` kept roughly the top 30% and killed the rest.

The cumulative sum is now walked to the number of cells that should *die*:

```python
n_kill = nb_tot - int(self.help_cutoff * nb_tot)
energy_idx = np.argmax(np.cumsum(sorted_nbs) >= n_kill)
```

Two further errors in the tie handling are fixed at the same time:

- `fraction_die_tie` divided the outstanding deaths by the *whole* surviving population rather
  than by the clones sitting exactly at the threshold, so the tie clones were scaled by a
  fraction that had nothing to do with their own size.
- The method pruned and recomputed `cut_idx` between the two stages, so the remainder was
  measured against an already-culled population and the tie clones took a second cull of
  roughly `1 - help_cutoff` on top of the first. The `JGF 12/14/23` comment in the file was
  someone noticing this; it is removed along with the intermediate `prune()`.

Checked against a pure-Python replica of both versions, ten clones of 100 cells with ΔG from
−20 to −11:

| `help_cutoff` | kept before | kept after |
| --- | --- | --- |
| 0.30 | 80% | 30% |
| 0.70 | 40% | 70% |
| 0.90 | 20% | 90% |

In every case the survivors are now strictly the tighter binders. Also guarded
`help_cutoff >= 1.0` (nothing dies) and an empty tie set, and replaced the seven per-cycle
debug prints with one summary line reporting cells kept, the fraction, and the ΔG cutoff.

> Note: this affects selection strength in **every** cycle of **every** run, so results
> generated on this branch before this fix — including the `codon-project/` CLK19gl runs — are
> not comparable to runs made after it. The pre-fix code is present on every feature branch and
> on `origin/main`, though `main` has it wrong in the opposite direction: `main` keeps 70% of
> cells but selects the *weakest* binders, since its kill mask is `binding_energies <
> energy_cutoff`.

## Known issues, not addressed

- `Mutability.csv` and `Substitution.csv` are re-read from disk and scanned linearly for every
  5-mer of every clone. Correct but slow; loading them once into a dict keyed by 5-mer would
  speed up large sweeps.
- `pick_memCells` uses `np.int`, removed in NumPy >= 1.24. Reached only on multi-injection runs.
- **`revert: true` is broken.** `run_recycle` calls `memory_cells.clone(residues, nbs)` at
  `lymph_node.py:553`, but `clone()` takes `(residues, nbs, parent_ids)` — a `TypeError` at the
  first recycle. Keep `revert: false`, the recommended setting, which is unaffected.
- The `memory_cells` placeholder built in `GerminalCenter.__init__` (`residues=None`, `nbs=1`,
  `E=inf`) is concatenated into the real memory pool on the first `append_mem_cells` call, so a
  phantom clone sits in the pool. Harmless today — nothing reads the pool on single-injection
  runs — but it must be cleared before prime-boost is revived, or `pick_memCells` can draw it.
- `config/__init__.py` uses `collections.Mapping`, removed in Python 3.10.
