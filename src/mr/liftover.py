"""Lift a GRCh38 summary-statistics table to GRCh37 (UCSC hg38ToHg19 chain).

    uv run --frozen python -m src.mr.liftover IN OUT --chrom COL --pos COL \\
        [--chain data/raw/mr/liftover/hg38ToHg19.over.chain.gz] [--log OUT.liftover_log.tsv]

pyliftover==0.4.1 is in uv.lock (added 2026-10-21). Standard library otherwise.

Output = every input column unchanged plus chrom_grch37 and pos_grch37
(1-based). A row is DROPPED - and counted by reason in the log - when the
position does not map, maps to several places, maps to the opposite strand
(alleles would need complementing; dropping is the conservative choice), maps
to a different chromosome, or has a non-numeric position. Nothing else is
changed.

Guard: a file under data/raw/mr/quarantine_outcome/ (Alzheimer's outcomes) is
refused unless the Module 1 pre-registration is frozen
(`python -m src.common.prereg check mr` succeeds), because lifting reads
every row.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
QUARANTINE = REPO / "data" / "raw" / "mr" / "quarantine_outcome"
DEFAULT_CHAIN = REPO / "data" / "raw" / "mr" / "liftover" / "hg38ToHg19.over.chain.gz"


def _open(path: Path, mode: str):
    return gzip.open(path, mode + "t", encoding="utf-8") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def _guard(path: Path) -> None:
    if QUARANTINE in path.resolve().parents:
        r = subprocess.run(["uv", "--directory", str(REPO), "run", "--frozen", "python", "-m", "src.common.prereg",
                            "check", "mr"], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"refused: {path} is a quarantined AD outcome and prereg/mr.md is not frozen.\n{r.stdout}{r.stderr}")


def lift(inp: Path, out: Path, chrom_col: str, pos_col: str, chain: Path, log: Path) -> Counter:
    from pyliftover import LiftOver  # imported here so --help works without it

    _guard(inp)
    lo = LiftOver(str(chain))
    n = Counter()
    with _open(inp, "r") as fi, _open(out, "w") as fo:
        header_line = fi.readline().rstrip("\n")
        sep = "\t" if "\t" in header_line else None
        header = header_line.split(sep)
        ic, ip = header.index(chrom_col), header.index(pos_col)
        osep = "\t"
        fo.write(osep.join(header + ["chrom_grch37", "pos_grch37"]) + "\n")
        for line in fi:
            n["read"] += 1
            f = line.rstrip("\n").split(sep)
            chrom = f[ic].replace("chr", "")
            chrom = {"23": "X", "24": "Y", "25": "M", "MT": "M"}.get(chrom, chrom)
            try:
                pos = int(float(f[ip]))
            except ValueError:
                n["dropped: position not numeric"] += 1
                continue
            hits = lo.convert_coordinate("chr" + chrom, pos - 1)  # pyliftover is 0-based
            if not hits:
                n["dropped: unmapped in chain"] += 1
                continue
            if len(hits) > 1:
                n["dropped: maps to several GRCh37 positions"] += 1
                continue
            c37, p37, strand, _ = hits[0]
            c37 = c37.replace("chr", "")
            if strand != "+":
                n["dropped: maps to the opposite strand"] += 1
                continue
            if c37 != chrom:
                n["dropped: maps to a different chromosome"] += 1
                continue
            fo.write(osep.join(f + [c37, str(p37 + 1)]) + "\n")
            n["written"] += 1
    with open(log, "w", encoding="utf-8") as fh:
        fh.write("input\toutput\tchain\toutcome\tn\n")
        for k, v in n.items():
            fh.write(f"{inp}\t{out}\t{chain}\t{k}\t{v}\n")
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--chrom", required=True, help="chromosome column (GRCh38)")
    ap.add_argument("--pos", required=True, help="1-based position column (GRCh38)")
    ap.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    ap.add_argument("--log", type=Path, default=None)
    a = ap.parse_args()
    log = a.log or a.output.with_name(a.output.name + ".liftover_log.tsv")
    n = lift(a.input, a.output, a.chrom, a.pos, a.chain, log)
    for k, v in n.items():
        print(f"{k}\t{v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
