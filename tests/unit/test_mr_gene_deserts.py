"""Unit tests for the Module 1 helpers added on 2026-10-21: the gene-desert
interval logic (src/mr/gene_deserts.py) and the standard-library .xlsx reader
(src/mr/xlsx.py). Hand-built fixtures only."""
import zipfile

from src.mr.gene_deserts import deserts, intergenic, merge
from src.mr.xlsx import iter_rows, sheet_names


def test_overlapping_and_touching_genes_are_merged():
    assert merge([(10, 20), (15, 30), (31, 40), (50, 60)]) == [(10, 40), (50, 60)]


def test_intervals_lie_between_two_genes_and_never_across_a_boundary_gap():
    genes = [(100, 200), (500, 600), (1000, 1100), (2000, 2100)]
    bounds = [(1, 1), (700, 800), (5000, 5000)]          # chromosome ends and one large gap
    ivs = intergenic(genes, bounds)
    # 201-499 kept; 601-999 crosses the gap -> excluded; 1101-1999 kept;
    # nothing before the first or after the last gene (edge stretches).
    assert [(lo, hi) for lo, hi, _, _ in ivs] == [(201, 499), (1101, 1999)]


def test_deserts_are_the_longest_fraction_with_ties_kept():
    ivs = {"1": [(1, 100, 0, 101), (201, 400, 200, 401), (501, 700, 500, 701)], "2": [(1, 50, 0, 51)]}
    chosen, n, cut = deserts(ivs, 0.25)                    # ceil(0.25 * 4) = 1 -> cut-off 200 bp, two tied
    assert n == 4 and cut == 200
    assert [(c, lo) for _, c, lo, _ in chosen] == [("1", 201), ("1", 501)]


def test_xlsx_reader_reads_shared_and_inline_strings_and_fills_gaps(tmp_path):
    p = tmp_path / "t.xlsx"
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rel = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook {ns} {rel}><sheets><sheet name="ST9" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", f'<sst {ns}><si><t>rsID</t></si><si><t>rs1611115</t></si></sst>')
        z.writestr("xl/worksheets/sheet1.xml",
                   f'<worksheet {ns}><sheetData>'
                   '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>BETA</t></is></c></row>'
                   '<row r="2"><c r="A2" t="s"><v>1</v></c><c r="C2"><v>1.1701900000000001</v></c></row>'
                   '</sheetData></worksheet>')
    assert sheet_names(p) == ["ST9"]
    assert list(iter_rows(p, "ST9")) == [["rsID", "", "BETA"], ["rs1611115", "", "1.1701900000000001"]]
