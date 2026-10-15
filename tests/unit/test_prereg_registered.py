"""The registered OSF documents kept in prereg/registered/ are the ones OSF serves.

Each PDF's sha256 must equal the hash OSF reported for the attached file
(prereg/registered/osf_registrations.yaml), and the transcription it names must exist.
"""
import pytest
import yaml

from src.common.paths import REPO_ROOT
from src.common.provenance import sha256_file

REG_DIR = REPO_ROOT / "prereg" / "registered"
REGS = yaml.safe_load((REG_DIR / "osf_registrations.yaml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("name", sorted(REGS))
def test_registered_pdf_is_byte_identical_to_osf(name):
    rec = REGS[name]
    pdf = REG_DIR / rec["file"]
    assert pdf.stat().st_size == rec["bytes"]
    assert sha256_file(pdf) == rec["sha256"]


@pytest.mark.parametrize("name", sorted(REGS))
def test_registration_is_public_not_withdrawn_and_transcribed(name):
    rec = REGS[name]
    assert rec["public"] is True and rec["withdrawn"] is False
    assert (REPO_ROOT / rec["transcription"]).is_file()
