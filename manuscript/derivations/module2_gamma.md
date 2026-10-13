# Module 2 — expected value of the amyloid effect γ (stated before any fit)

> **Status: derived by Claude; MACHINE-CHECKED 2026-10-13, NOT human-checked.**
> `checks/module2_gamma.py` checks the algebra symbolically (including the
> age-adjusted form), the numbers, and claim A3 numerically. For A3, an
> individual-level model with HR 2.43 (5.0) is fitted by the population form
> as HR 2.13 (3.26): attenuated towards 1, as stated. 9/9 checks pass.
>
> **Role since the 2026-10-13 revision: descriptive reference only.** Module 2
> no longer tests amyloid gating (prereg/baserate.md change log). The fitted
> HR_A comes from the *age-adjusted* variant A,
> μ(t) = m(t)(1 + (HR_A − 1)A(t)) with m(t) = exp(α + β(t − 50)/10). The
> mapping HR_A = amyloid-positive / amyloid-negative rate ratio at the same
> age is unchanged (§1 with μ0 replaced by m(t)). §5's comparison is
> reported only as a description, never as a verdict.
>
> The result still depends heavily on an untested **proxy assumption** (§3,
> A1).

## 1. What γ is in the model

Variant A gives the LC → cortex transition (stage c−1 → c, i.e. subcortical
stage c → cortical stage 1a) the population rate

    μ(t) = μ0 + γ A(t),

where A(t) is the probability that a person of age t is amyloid-positive.
Interpret it as a two-type mixture: an amyloid-negative person makes the
transition at rate μ0, an amyloid-positive person at rate μ0 + γ, and a
fraction A(t) of people are positive:

    μ(t) = (1 − A) μ0 + A (μ0 + γ) = μ0 + γ A.

The amyloid-positive vs amyloid-negative **hazard ratio** for this
transition is therefore

    HR_A = (μ0 + γ) / μ0 = 1 + γ/μ0,   so   γ = μ0 (HR_A − 1).

The model is parameterised as μ0 = exp(log_mu0) and HR_A = exp(log_hr_amyloid)
(src/baserate/model.py). This is the same model as the brief's
μ(t) = μ0 + γA(t); it keeps every rate positive for any parameter value
(HR_A < 1 is allowed, so the fit can contradict the hypothesis) and makes
the expected value directly comparable with a published hazard ratio.

γ itself is in per-year units and scales with μ0, which is only known after
the fit. The expectation is therefore stated for the dimensionless HR_A
(equivalently γ/μ0); the implied γ is reported after the fit as
μ0 × (HR_A − 1).

## 2. Published hazard ratio used

Ossenkoppele R, Pichet Binette A, Groot C, et al. *Nat Med* 2022;28:2381–2387
(doi:10.1038/s41591-022-02049-x), **Table 2**, all ages: cognitively
unimpaired, tau-PET-negative participants, amyloid-PET-positive (A+T−,
n = 292) vs amyloid-negative (A−T−, n = 781), Cox model for progression to
MCI adjusted for age, sex, education and cohort:

    HR = 2.43, 95% CI 1.38–4.26.

(Extracted to data/extracted/amyloid_hr_ossenkoppele2022.csv,
`verified: false`.)

## 3. Mapping assumptions (each one could be wrong)

**A1 (proxy, the big one).** The amyloid HR for *clinical progression to MCI
among tau-PET-negative people* equals the amyloid HR for the *LC → cortex
tau transition*. Justification: the A+T− vs A−T− contrast holds detectable
cortical tau fixed (negative) and varies amyloid, and progression from that
state to MCI requires cortical tau to spread. Against: tau-PET cannot see
the LC or Braak stages 1a–II, so "T−" mixes people before and after the
modelled transition; and progression to MCI involves later transitions
too. I searched for an amyloid-stratified hazard ratio for *incident tau
positivity* (closer to the modelled transition) and found none in open
access (data/extracted/SOURCES.md).

**A2 (proportional hazards).** The HR does not vary with age. The paper's own
age strata (Table 2: 0.94 at 50–69, 2.93 at 70–79, 2.41 at 80+) are
compatible with this within their wide CIs.

**A3 (mixture at the rate level).** μ(t) = μ0 + γA(t) treats amyloid status
as independent of stage among people still at stage c. In reality
amyloid-positive people leave stage c faster, so those still there are
depleted of amyloid-positives, and the population A(t) overstates the
amyloid-positive share at risk. The fit then needs a smaller HR_A to explain
the same transitions: **the fitted HR_A is expected to be attenuated
towards 1** relative to an individual-level HR. A fitted HR_A below 2.43 is
therefore weak evidence against the expectation; a fitted HR_A ≤ 1 is not
explained by attenuation.

## 4. Expected value

On the log scale (the model's parameter):

    log HR_A = ln 2.43 = 0.888
    95% CI:   ln 1.38 = 0.322  to  ln 4.26 = 1.449
    implied SE = (1.449 − 0.322) / (2 × 1.96) = 0.288

Equivalently **HR_A = 2.43 (1.38–4.26)**, γ/μ0 = 1.43 (0.38–3.26), and
γ = 1.43 μ0 per year once μ0 is fitted.

## 5. Comparison rule (fixed now, before the fit)

Let I_pub = [0.322, 1.449] (published 95% CI for log HR) and I_fit = the
fit's 95% credible interval for log_hr_amyloid.

* **consistent**: 0.888 lies inside I_fit;
* **tension**: 0.888 lies outside I_fit but the intervals overlap;
* **disagrees** (flagged loudly in the run log and report): I_pub and I_fit
  do not overlap.

Because of A1 and A3, the report will describe any of these outcomes as
dependent on an unconfirmed proxy and a derivation the researcher has not
yet checked.
