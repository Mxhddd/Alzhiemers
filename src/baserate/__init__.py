"""Module 2 — Markov base-rate model."""

import os

# Let JAX expose one CPU device per core so NUTS chains run in parallel. Must be
# set before JAX is first imported (src/baserate imports JAX lazily).
os.environ.setdefault("XLA_FLAGS", f"--xla_force_host_platform_device_count={os.cpu_count() or 1}")
