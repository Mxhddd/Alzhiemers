# Pre-registration DRAFT — Module 6: network diffusion (LC vs entorhinal seed)

> **DRAFT by Claude for researcher review; not frozen.** Save as
> `prereg/network.md`, upload to OSF, `make freeze-prereg NAME=network OSF=<url>`.

## Prediction and what counts against it
The LC seed fits the earliest-stage regional tau pattern better than the
entorhinal seed: Δr = r_LC − r_EC > 0 with spin p < [0.05]. Against: Δr ≤ 0,
or the difference is not distinguishable from the spin null, or both seeds'
best fits lie beyond the crossover time (then the data cannot discriminate —
reported as uninformative, not as support).

## Methods fixed in advance
- Connectome [HCP dataset, parcellation]; tau [ADNI, earliest Braak-stage group, N].
- Normalised Laplacian; τ grid log-spaced 0.01–50 (400 points); each seed's r
  maximised over τ.
- Significance ONLY from spin nulls (10,000 rotations, Hungarian reassignment;
  right hemisphere mirrored), with the τ-maximisation repeated inside the null.
- Controls: [seeds]. Crossover threshold 0.95.
- The spin test was validated on a known null: 6% false positives (naive
  shuffling: 57%) — tests/validation/test_network.py.
