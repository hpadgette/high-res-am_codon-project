# Atomistic Affinity Maturation Model

This is a preliminary repository for the atomistic affinity maturation model (https://doi.org/10.1371/journal.pcbi.1009391). This will be undergoing some refactoring and reorganizing over the coming weeks, and shouldn't be used until this readme is updated. This readme will serve as a living document throughout these changes

# Plan:
- Add a configuration file (yaml/json dictionary format)
- Remove the python files which have any sort of sequence information (move to experiment json)
- Remove parameters that are set in the `mixture.py` file (move to common json)
- Add a task/pathogen flag and directory so the same code can be used for any disease
- Add sum of single point (SSP) functionality for binding affinity scoring
- Update `main` to use `ppdx` rather than `ppdg`
