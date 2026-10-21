"""Check prereg/<name>_osf_answers.txt (answers for the OSF Preregistration form).

    python -m src.registration.osf_answers <name>      (exit 1 on any problem)

Format of the answers file (plain text, no markdown):

    ATTACH prereg/<name>.pdf TO THESE FILE BOXES: Study design; Statistical models

    STEP 1 - METADATA

    Title:
    <answer, one or more lines>

    Description:
    <answer>
    ...

Every field label below must appear, in this order, on a line of its own
ending in a colon, followed by a non-empty answer. Checks:
  * all fields present, in order, each answered;
  * the attach line is the first line and names prereg/<name>.pdf;
  * Foreknowledge quotes exactly one of the OSF options, verbatim;
  * Study type / Intention for causal interpretation / Blinding use an OSF option;
  * License is CC0 1.0 Universal;
  * no markdown (no **, __, backticks, '#' headings, '|' tables);
  * every number in the answers also appears in prereg/<name>.md
    (ignoring numbers inside URLs, the license, and the STEP numbers);
  * the context field links both existing registrations and discloses Claude.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PREREG = REPO_ROOT / "prereg"

STEPS = {
    "STEP 1 - METADATA": ["Title", "Description", "License", "Subjects", "Tags"],
    "STEP 2 - OVERVIEW": ["Research questions or hypotheses", "Foreknowledge of data or evidence",
                          "Explanation of foreknowledge and managing unintended influences"],
    "STEP 3 - RESEARCH DESIGN": ["Study type", "Intention for causal interpretation",
                                 "Blinding of experimental treatments", "Additional blinding during research or analysis",
                                 "Study design", "Randomization"],
    "STEP 4 - SAMPLING": ["Data collection procedures", "Sample size", "Sample size rationale",
                          "Starting and stopping rules"],
    "STEP 5 - VARIABLES": ["Manipulated variables", "Measured variables", "Indices"],
    "STEP 6 - ANALYSIS PLAN": ["Statistical models", "Transformations", "Inference criteria",
                               "Data inclusion and exclusion", "Missing data", "Other planned analysis"],
    "STEP 7 - OTHER": ["Context and additional information"],
}
FOREKNOWLEDGE = [
    "Data does not yet exist.",
    "Data exists but the authors cannot observe it yet.",
    "Data exists but the authors have not observed it yet.",
    "Only people other than the authors have observed the data.",
    "Authors' limited observation of the data could not influence their analysis decisions.",
    "Authors have observed the data, but have not performed the proposed analyses.",
    "Authors have observed the data.",
    "Analyses in this plan have been conducted already.",
]
STUDY_TYPES = ["Randomized Experiment", "Non-randomized study", "Systematic review",
               "Adjustment or control-based design", "Quasi-experimental design", "Descriptive study",
               "Simulation study", "Other"]
CAUSAL = ["No causal relationship inferred", "Indirect inference", "Direct inference"]
EXISTING = ["https://osf.io/xb7fy/", "https://osf.io/9kb4w/"]
NUM = re.compile(r"(?<![\w.])[-−]?\d[\d,]*(?:\.\d+)?(?:[eE][-−]?\d+)?(?:\s?[×x]\s?10[⁻\-−]?[⁰¹²³⁴⁵⁶⁷⁸⁹\d]+)?")


def parse(text: str) -> tuple[str, dict[str, str], list[str]]:
    lines = text.splitlines()
    labels = [lab for step in STEPS.values() for lab in step]
    order, answers, cur, buf = [], {}, None, []
    for ln in lines[1:]:
        s = ln.strip()
        if s in STEPS:
            continue
        if s.endswith(":") and s[:-1] in labels:
            if cur:
                answers[cur] = "\n".join(buf).strip()
            cur, buf = s[:-1], []
            order.append(cur)
        elif cur:
            buf.append(ln)
    if cur:
        answers[cur] = "\n".join(buf).strip()
    return (lines[0] if lines else ""), answers, order


def _nums(s: str) -> set[str]:
    s = re.sub(r"https?://\S+", " ", s)
    s = s.replace("CC0 1.0", " ")
    out = set()
    for m in NUM.finditer(s):
        t = m.group(0).replace(",", "").replace("−", "-").lstrip("-").replace(" ", "")
        out.add(t)
    return out


def _doc_nums(doc: str) -> set[str]:
    base = _nums(doc)
    extra = set()
    for t in base:   # 1,000 written as 1000 and vice versa are the same number
        extra.add(t.rstrip("0").rstrip(".") if "." in t else t)
    return base | extra


def check(name: str) -> list[str]:
    path = PREREG / f"{name}_osf_answers.txt"
    md = PREREG / f"{name}.md"
    if not path.exists():
        return [f"{path} does not exist"]
    if not md.exists():
        return [f"{md} does not exist"]
    text = path.read_text(encoding="utf-8")
    problems = []
    first, answers, order = parse(text)
    want_first = f"ATTACH prereg/{name}.pdf TO THESE FILE BOXES: Study design; Statistical models"
    if first.strip() != want_first:
        problems.append(f"first line must be exactly: {want_first!r}")
    for step in STEPS:
        if step not in text:
            problems.append(f"missing step heading {step!r}")
    labels = [lab for step in STEPS.values() for lab in step]
    if order != labels:
        missing = [lab for lab in labels if lab not in order]
        problems.append(f"fields missing or out of order; missing: {missing}" if missing else "fields out of order")
    for lab in labels:
        if lab in answers and not answers[lab].strip():
            problems.append(f"field {lab!r} has no answer")
    fk = answers.get("Foreknowledge of data or evidence", "").strip()
    if fk not in FOREKNOWLEDGE:
        problems.append(f"Foreknowledge must be exactly one OSF option, verbatim; got {fk!r}")
    st = answers.get("Study type", "").strip().splitlines()[0] if answers.get("Study type") else ""
    if st.rstrip(".") not in STUDY_TYPES:
        problems.append(f"Study type must start with one OSF option on its own line; got {st!r}")
    ci = answers.get("Intention for causal interpretation", "").strip().splitlines()[0] \
        if answers.get("Intention for causal interpretation") else ""
    if ci.rstrip(".") not in CAUSAL:
        problems.append(f"Intention for causal interpretation must start with one OSF option; got {ci!r}")
    if answers.get("License", "").strip() != "CC0 1.0 Universal":
        problems.append("License must be exactly 'CC0 1.0 Universal'")
    if re.search(r"\*\*|__|`|^\s*#|^\s*\|", text, flags=re.M):
        problems.append("markdown formatting found (**, __, backticks, # headings or | tables)")
    ctx = answers.get("Context and additional information", "")
    for url in EXISTING:
        if url not in ctx:
            problems.append(f"Context field must link {url}")
    if "Claude" not in ctx or "Anthropic" not in ctx:
        problems.append("Context field must disclose the use of Claude (Anthropic)")
    body = "\n".join(answers.values())
    doc_nums = _doc_nums(md.read_text(encoding="utf-8"))
    stray = sorted(n for n in _nums(body) if n not in doc_nums and
                   (n.rstrip("0").rstrip(".") if "." in n else n) not in doc_nums)
    if stray:
        problems.append(f"numbers in the answers that do not appear in prereg/{name}.md: {stray}")
    if re.search(r"@[A-Za-z0-9-]+\.[A-Za-z]{2,}", text):
        problems.append("an email address appears in the answers file")
    return problems


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m src.registration.osf_answers <name>")
    probs = check(sys.argv[1])
    for p in probs:
        print("PROBLEM:", p)
    print("OK" if not probs else f"{len(probs)} problem(s)")
    sys.exit(1 if probs else 0)
