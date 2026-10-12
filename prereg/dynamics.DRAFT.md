# Pre-registration DRAFT — Module 4: ODE dynamical model

> **DRAFT by Claude for researcher review; not frozen.** Fill brackets, save
> as `prereg/dynamics.md`, upload to OSF, `make freeze-prereg NAME=dynamics OSF=<url>`.

## Questions
1. Can one parameter set make both atomoxetine and prazosin protective?
2. Is the system bistable in the second-hit parameter [name]?
3. Does the dose-response show a U-shape?

## Model
Equations: [researcher's final equations; see manuscript/notes/module4_model_proposal.md].
Parameters: config/dynamics_params.yaml at commit [hash]; weak parameters swept
over their stated ranges, never reported as fitted values.

## Pre-specified outputs and interpretation
- U-shape: both dose-range ends exceed the interior minimum by > 5%.
- Robustness: fraction of the Latin-hypercube sample (n = 1000) with a U-shape
  and with both drugs protective, with Wilson 95% CIs. **These fractions are
  reported whatever they are.** [Researcher: state in advance what fraction
  would count as "robust" vs "fragile", e.g. < 10% = fragile.]
- Bistability: >= 2 stable fixed points over some interval of the scan.
- A model in which the answers depend on solver or tolerance (Radau vs BDF,
  tolerances tightened 10x, disagreement > 0.1%) is not interpreted.
