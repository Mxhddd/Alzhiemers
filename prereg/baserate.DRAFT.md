# Pre-registration DRAFT — Module 2: Markov base-rate model

> **Status: DRAFT written by Claude from PROJECT_BRIEF.md for the researcher to
> review. Not frozen.** Fill the bracketed items, save as `prereg/baserate.md`,
> upload to OSF, `make freeze-prereg NAME=baserate OSF=<url>`, commit.

## Question

If nearly everyone has LC tau by 40, what must be true for only a minority to
develop Alzheimer's dementia?

## Models (all fitted; none favoured)

- **A, gated:** LC → cortex rate μ(t) = μ₀ + γ·A(t), A(t) = amyloid-positivity
  probability (piecewise constant by age band, researcher-supplied).
- **B, independent:** cortical tau arises at rate ν(t) from any pre-cortical
  stage, regardless of LC stage.
- **C, multi-hit:** cortical entry after m independent hits (m = [2]).

Other transitions: q(t) = exp(α + β(t − 50)).

## Data

Braak stage counts by age band [sources], transcribed by the researcher.
Stage order [ ]; first cortical stage [ ]; dementia threshold stage [ ].

## Predictions and what would count against the hypothesis

- Supportive: A has the best cross-validated held-out log-likelihood, and γ > 0
  with its 95% credible interval excluding 0.
- Against: γ's interval includes 0 or is negative, or B/C fit as well or
  better by cross-validation [define margin, e.g. Δ held-out log-lik > 2].
- Falsification (any model): predicted dementia incidence differs from the
  published incidence by more than [tolerance on |log ratio|] in any age band.
- Researcher's hand-derived γ: [value, source]; disagreement is flagged.

## Analysis settings

Matrix-exponential step [h years, from the hand error analysis]; MLE with 20
starts; NUTS 4 × 1000 draws; R-hat ≤ 1.01 and zero divergences required for
the posterior to be reported without caveat; 5-fold CV by age band.
