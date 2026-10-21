#!/usr/bin/env bash
# Keep the registered Module 2 run going across container restarts.
# Prints one of: RUNNING, COMPLETE, FAILED (with the error), RELAUNCHED.
# Relaunching is safe: the run resumes from its stage and sub-stage
# checkpoints, which are keyed to config + code + settings + inputs.
cd "$(dirname "$0")/.." || exit 2
LOGDIR=results/logs/real
latest=$(ls -t $LOGDIR/baserate_run*.log 2>/dev/null | head -1)
if ps -eo cmd | grep -v "bash -c" | grep -q "[p]ython3 -m src.baserate.run_baserate"; then
  echo "RUNNING ($latest)"; exit 0
fi
if [ -n "$latest" ] && grep -q "ALL STAGES COMPLETE" "$latest"; then
  echo "COMPLETE ($latest)"; exit 0
fi
if [ -n "$latest" ] && grep -qE "Traceback|Error:|RuntimeError" "$latest"; then
  echo "FAILED ($latest):"; grep -E "Traceback|Error" "$latest" | tail -3; exit 1
fi
n=$(ls $LOGDIR/baserate_run*.log 2>/dev/null | wc -l); new="$LOGDIR/baserate_run$((n+1)).log"
setsid nohup uv run --frozen python -m src.baserate.run_baserate > "$new" 2>&1 < /dev/null & disown
echo "RELAUNCHED -> $new"
