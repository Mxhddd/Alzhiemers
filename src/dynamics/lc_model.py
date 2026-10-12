"""The LC rate-controller ODE model — RESEARCHER TO SPECIFY.

The brief assigns the mathematics to the researcher. This file defines the
interface the machinery needs; every function raises until the researcher
writes it. A candidate structure (equations without numbers) is proposed in
manuscript/notes/module4_model_proposal.md for the researcher to accept, change or reject.

Interface:
  build_model() -> ODEModel            the right-hand side, state names, non-negative
                                       states and conserved quantities
  drugs() -> dict[str, Drug]           at least "atomoxetine" and "prazosin": how a
                                       dose changes parameters
  outcome(traj) -> float               the scalar the dose-response is judged on
                                       (lower = better), e.g. cortical tau at final time
  initial_state(p) -> np.ndarray       starting state
  coupling_params() -> list[str]       parameters that couple subsystems; setting them
                                       to zero must reduce the model to known exponentials
  uncoupled_solution(t, p, y0) -> np.ndarray
                                       the researcher's closed-form solution of that
                                       degenerate case (validation, brief §5)
Parameters live in config/dynamics_params.yaml with units, source and confidence.
"""

from __future__ import annotations

_MSG = ("The LC model has not been specified. The researcher writes the equations "
        "(see manuscript/notes/module4_model_proposal.md) and parameters (config/dynamics_params.yaml).")


def build_model():
    raise NotImplementedError(_MSG)


def drugs():
    raise NotImplementedError(_MSG)


def outcome(traj):
    raise NotImplementedError(_MSG)


def initial_state(p):
    raise NotImplementedError(_MSG)


def coupling_params():
    raise NotImplementedError(_MSG)


def uncoupled_solution(t, p, y0):
    raise NotImplementedError(_MSG)
