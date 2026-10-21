"""Every registration package in prereg/ passes its consistency checks."""
import pytest

from src.registration import osf_answers
from src.registration.pdf import PREREG, text_layer_ok

NAMES = sorted(p.name[:-len("_osf_answers.txt")] for p in PREREG.glob("*_osf_answers.txt"))


@pytest.mark.parametrize("name", NAMES)
def test_osf_answers_are_complete_and_consistent_with_the_document(name):
    assert osf_answers.check(name) == []


@pytest.mark.parametrize("name", NAMES)
def test_pdf_exists_and_has_a_text_layer_matching_the_document(name):
    pdf = PREREG / f"{name}.pdf"
    assert pdf.exists()
    assert pdf.stat().st_mtime >= (PREREG / f"{name}.md").stat().st_mtime - 1, "PDF older than the document: rebuild"
    ok, missing = text_layer_ok(PREREG / f"{name}.md", pdf)
    assert ok, missing


@pytest.mark.parametrize("name", NAMES)
def test_document_has_no_draft_status_block_and_names_the_author(name):
    text = (PREREG / f"{name}.md").read_text(encoding="utf-8")
    head = text[:1500].lower()
    assert "mohammad hassan" in head
    # no status block (the approved Tools wording "draft this text" is fine)
    for phrase in ("status:", "complete draft", "not yet frozen", "draft — ", "draft, ready", "do not freeze"):
        assert phrase not in head, phrase
    assert "## tools" in text.lower() and "claude" in text.lower() and "change log" in text.lower()


def test_checker_catches_problems(tmp_path, monkeypatch):
    monkeypatch.setattr(osf_answers, "PREREG", tmp_path)
    (tmp_path / "x.md").write_text("# T\n\nWe use 5 folds and a margin of 2.\n", encoding="utf-8")
    good = ["ATTACH prereg/x.pdf TO THESE FILE BOXES: Study design; Statistical models", ""]
    for step, labs in osf_answers.STEPS.items():
        good += [step, ""]
        for lab in labs:
            ans = {"License": "CC0 1.0 Universal",
                   "Foreknowledge of data or evidence": "Data does not yet exist.",
                   "Study type": "Simulation study", "Intention for causal interpretation": "No causal relationship inferred",
                   "Context and additional information":
                       "Claude (Anthropic) helped. https://osf.io/xb7fy/ https://osf.io/9kb4w/"}.get(lab, "5 folds, margin 2.")
            good += [f"{lab}:", ans, ""]
    (tmp_path / "x_osf_answers.txt").write_text("\n".join(good), encoding="utf-8")
    assert osf_answers.check("x") == []
    bad = "\n".join(good).replace("5 folds, margin 2.", "7 folds, **margin** 2.", 1)
    (tmp_path / "x_osf_answers.txt").write_text(bad, encoding="utf-8")
    probs = osf_answers.check("x")
    assert any("markdown" in p for p in probs) and any("'7'" in p for p in probs)
