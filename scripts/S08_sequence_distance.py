#!/usr/bin/env python3
"""Rebuild Supplementary Figure S8 from retained local sequence sources.

Requires Pillow and matplotlib (for portable font lookup). Run with:
    python3 scripts/S08_sequence_distance.py

The DUSA row is extracted from the saved CP070274.1 GenBank record, rather
than being copied from OTU2. The saved extracted reference FASTA and the
existing forced 16S alignment are independently checked before rendering.
"""

import csv
import hashlib
import itertools
import json
import re
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from matplotlib.font_manager import findfont


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/Tim16S.OTu.fas"
ALIGNMENT = ROOT / "output/phylogeny/argiope_16S_457_forced_backbone.fasta"
REFERENCE = ROOT / "output/species_symbionts/reference_validation/CP070274.genbank"
REFERENCE_FASTA = REFERENCE.with_name("CP070274_16S.fasta")
DATA = ROOT / "output/symbiont_16S_distance"
CSV = ROOT / "output/tables/S08_symbiont_16S_distance_data.csv"
FIGURE = ROOT / "figures/supplement/S08_symbiont_16S_distance.png"
IDS = ("CP070274.1_16S", "OTU2", "OTU287", "OTU3", "OTU145")
LABELS = {
    "CP070274.1_16S": "DUSA reference",
    "OTU2": "DUSA-like OTU2",
    "OTU287": "Close sequence relative OTU287",
    "OTU3": "Mesoplasma OTU3",
    "OTU145": "Mesoplasma OTU145",
}
NAVY = "#19394b"
MUTED = "#546c79"
GRID = "#d8e2e8"
PALE = "#f4f7f8"
BASE_BG = {"A": "#d8eee6", "C": "#dceafa", "G": "#fce9ce", "T": "#f4dde9"}
BASE_FG = {"A": "#13674f", "C": "#245d98", "G": "#995908", "T": "#93416b"}
REGULAR = os.environ.get("ARGIOPE_SANS_FONT", findfont("DejaVu Sans"))
MONO = os.environ.get("ARGIOPE_MONO_FONT", findfont("DejaVu Sans Mono"))


def read_fasta(path):
    result = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0].upper()
            assert name not in result, f"Duplicate FASTA identifier: {name}"
            result[name] = ""
        elif name:
            result[name] += line.strip().upper()
    return result


def deposited_16s():
    """Read the single simple complementary 16S feature in this local record.

    This deliberately rejects a changed feature structure; no general GenBank
    parsing or guessed extraction coordinates are used.
    """
    record = REFERENCE.read_text()
    accession = re.search(r"^VERSION\s+(\S+)", record, re.MULTILINE).group(1)
    assert accession == "CP070274.1"
    features = record.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0]
    blocks = re.split(r"(?m)(?=^     \S)", features)
    matches = [b for b in blocks if b.startswith("     rRNA ")
               and '/product="16S ribosomal RNA"' in b]
    assert len(matches) == 1, "Expected a single annotated 16S rRNA feature"
    feature = re.search(r"^     rRNA\s+complement\((\d+)\.\.(\d+)\)", matches[0])
    assert feature, "Unsupported GenBank 16S feature structure"
    start, end = map(int, feature.groups())
    genome = "".join(re.findall(r"[acgt]+", record.split("\nORIGIN", 1)[1].split("//", 1)[0])).upper()
    rrna = genome[start - 1:end].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert len(rrna) == end - start + 1 == 1530
    saved = read_fasta(REFERENCE_FASTA)
    assert rrna == saved["CP070274.1_16S"], "Saved reference FASTA disagrees with GenBank feature"
    return accession, rrna, start, end


def sequences_and_provenance():
    raw = read_fasta(RAW)
    aligned = read_fasta(ALIGNMENT)
    otus = IDS[1:]
    assert len({len(aligned[key]) for key in otus}) == 1
    # The four local fragments occupy exactly the same 101 backbone columns.
    # Removing only their shared gaps retains homologous sites for every pair.
    retained = [i for i in range(len(aligned["OTU2"]))
                if any(aligned[key][i] != "-" for key in otus)]
    assert len(retained) == 101
    local = {key: "".join(aligned[key][i] for i in retained) for key in otus}
    assert all(len(seq) == 101 and set(seq) <= set("ACGT") for seq in local.values())
    assert all(local[key] == raw[key] for key in otus), "Local FASTA and aligned fragments disagree"
    accession, rrna, gene_start, gene_end = deposited_16s()
    offset = rrna.find(local["OTU2"])
    assert offset >= 0 and rrna.find(local["OTU2"], offset + 1) < 0
    # The reference bases come from the deposited gene, anchored by the uniquely
    # matching local fragment. No bases outside the measured overlap are added.
    reference_fragment = rrna[offset:offset + 101]
    sequences = {IDS[0]: reference_fragment, **local}
    assert list(sequences) == list(IDS)
    paths = [RAW, ALIGNMENT, REFERENCE, REFERENCE_FASTA]
    provenance = {
        "figure": str(FIGURE.relative_to(ROOT)),
        "script": "scripts/S08_sequence_distance.py",
        "reference_accession": accession,
        "reference_gene": "16S ribosomal RNA",
        "reference_gene_genome_start_1based": gene_start,
        "reference_gene_genome_end_1based": gene_end,
        "reference_gene_strand": "minus",
        "reference_gene_length_bp": len(rrna),
        "reference_overlap_start_1based_in_oriented_16S": offset + 1,
        "reference_overlap_end_1based_in_oriented_16S": offset + 101,
        "aligned_sites_bp": len(retained),
        "retained_backbone_columns_1based": [i + 1 for i in retained],
        "alignment_method": "Retain the 101 common ACGT sites of the existing forced backbone alignment; deposited DUSA fragment anchored by unique exact OTU2 match",
        "distance_method": "Uncorrected p-distance: mismatches / comparable ACGT sites",
        "ambiguity_or_gap_sites_in_retained_columns": 0,
        "local_OTU_labels": {key: LABELS[key] for key in otus},
        "limitation": "OTU2 matches the deposited DUSA reference across this 101-bp fragment only. OTU287 is a distinct close sequence relative (4/101 differences). These fragment distances do not establish identical strains, whole genomes, infection status, or species boundaries; no species cutoff is asserted.",
        "sources": {str(p.relative_to(ROOT)): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in paths},
    }
    return sequences, provenance


def distances(sequences, provenance):
    result = []
    for a, b in itertools.combinations(IDS, 2):
        pairs = [(x, y) for x, y in zip(sequences[a], sequences[b]) if x in "ACGT" and y in "ACGT"]
        n = len(pairs)
        differences = sum(x != y for x, y in pairs)
        assert n == 101
        result.append({
            "sequence_a": a, "sequence_b": b,
            "label_a": LABELS[a], "label_b": LABELS[b],
            "comparable_aligned_ACGT_sites": n,
            "mismatches": differences,
            "uncorrected_p_distance": differences / n,
            "distance_percent": 100 * differences / n,
            "identity_percent": 100 * (1 - differences / n),
            "reference_accession": provenance["reference_accession"],
            "reference_overlap_start_1based_in_oriented_16S": provenance["reference_overlap_start_1based_in_oriented_16S"],
            "reference_overlap_end_1based_in_oriented_16S": provenance["reference_overlap_end_1based_in_oriented_16S"],
        })
    return result


def compact_columns(sequences):
    variable = {i for i in range(101) if len({sequences[key][i] for key in IDS}) > 1}
    columns = []
    i = 0
    while i < 101:
        if i in variable:
            columns.append((i, i, False))
            i += 1
            continue
        end = i
        while end + 1 < 101 and end + 1 not in variable:
            end += 1
        if end - i + 1 >= 3:
            columns.append((i, end, True))
        else:
            columns.extend((j, j, False) for j in range(i, end + 1))
        i = end + 1
    assert [j for start, end, _ in columns for j in range(start, end + 1)] == list(range(101))
    return columns


class Canvas:
    def __init__(self, width, height):
        self.image = Image.new("RGB", (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.fonts = {}

    def text(self, x, y, value, size, fill=NAVY, mono=False, anchor=None):
        key = (size, mono)
        if key not in self.fonts:
            self.fonts[key] = ImageFont.truetype(MONO if mono else REGULAR, size)
        self.draw.text((x, y), value, font=self.fonts[key], fill=fill, anchor=anchor)

    def box(self, xy, fill, outline=None, radius=0, width=2):
        self.draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)

    def line(self, xy, fill, width=2):
        self.draw.line(xy, fill=fill, width=width)


def render(sequences, rows, provenance):
    width, height = 4000, 2440
    d = Canvas(width, height)
    d.text(72, 55, "A   16S sequence alignment", 54)
    reference_start = provenance["reference_overlap_start_1based_in_oriented_16S"]
    reference_end = provenance["reference_overlap_end_1based_in_oriented_16S"]
    d.text(72, 132, f"101 homologous sites; DUSA reference positions {reference_start}–{reference_end}", 36, MUTED)
    top, row_step = 272, 123
    x = 640
    for start, end, is_range in compact_columns(sequences):
        cw = 127 if is_range else 66
        center = x + cw / 2
        coord = str(start + 1) if start == end else f"{start + 1}–{end + 1}"
        d.text(center, top - 53, coord, 30, MUTED, mono=True, anchor="mm")
        for row, key in enumerate(IDS):
            y = top + row * row_step
            if is_range:
                d.box((x + 3, y, x + cw - 3, y + 75), "#e7ecee", radius=7)
                d.text(center, y + 37, "…", 40, MUTED, mono=True, anchor="mm")
            else:
                base = sequences[key][start]
                variable = len({sequences[k][start] for k in IDS}) > 1
                d.box((x + 3, y, x + cw - 3, y + 75), BASE_BG[base] if variable else PALE, GRID, radius=7)
                d.text(center, y + 37, base, 43, BASE_FG[base] if variable else MUTED, mono=True, anchor="mm")
        x += cw
    assert x < width - 70, (x, width)
    d.text(72, top - 5, LABELS[IDS[0]], 43)
    d.text(72, top + 44, "CP070274.1", 33, MUTED)
    for row, key in enumerate(IDS[1:], 1):
        if key == "OTU287":
            d.text(72, top + row * row_step - 5, "Close sequence relative", 41)
            d.text(72, top + row * row_step + 44, "OTU287", 33, MUTED)
        else:
            d.text(72, top + row * row_step + 13, LABELS[key], 41)
    d.text(72, 901, "Colored bases vary among sequences; gray ranges collapse ≥3 identical bases. Positions are relative to the fragment.", 33, MUTED)
    d.line((72, 969, width - 72, 969), GRID, 3)
    d.text(72, 1016, "B   Pairwise sequence distance", 54)
    d.text(72, 1094, "Uncorrected differences / 101 aligned ACGT sites", 36, MUTED)
    chart_left, chart_width, chart_top = 2100, 1150, 1251
    chart_bottom = chart_top + (len(rows) - 1) * 99 + 48
    for tick in (0, 5, 10, 15, 20):
        tx = chart_left + chart_width * tick / 20
        d.line((tx, chart_top - 42, tx, chart_bottom), "#e1e9ed", 2)
        d.text(tx, chart_top - 80, f"{tick}%", 31, MUTED, anchor="mt")
    for row, data in enumerate(rows):
        y = chart_top + row * 99
        d.text(72, y - 7, f"{data['label_a']}  /  {data['label_b']}", 41)
        end = chart_left + chart_width * data["distance_percent"] / 20
        if data["mismatches"]:
            color = "#62958b" if data["sequence_a"] == "OTU3" else "#5d88a6"
            d.box((chart_left, y - 3, end, y + 40), color, radius=7)
        else:
            d.draw.ellipse((chart_left - 6, y + 12, chart_left + 6, y + 24), fill=NAVY)
        d.text(3340, y - 7, f"{data['mismatches']}/101  ({data['distance_percent']:.2f}%)", 41)
    d.line((72, 2233, width - 72, 2233), GRID, 3)
    d.text(72, 2261, "DUSA reference and OTU2 match in this 101-bp fragment; OTU287 is a distinct close relative.", 36, MUTED)
    d.text(72, 2315, "Fragment distances do not resolve strain identity or species boundaries.", 36, MUTED)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    d.image.save(FIGURE, dpi=(300, 300), optimize=True)


def write_outputs(sequences, rows, provenance):
    DATA.mkdir(parents=True, exist_ok=True)
    CSV.parent.mkdir(parents=True, exist_ok=True)
    with CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (DATA / "aligned_sequences.fasta").write_text("".join(
        f">{key} {LABELS[key]}; 101 homologous fragment sites\n{sequences[key]}\n" for key in IDS))
    (DATA / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    reference_start = provenance["reference_overlap_start_1based_in_oriented_16S"]
    reference_end = provenance["reference_overlap_end_1based_in_oriented_16S"]
    lines = ["# Supplementary Figure S8: symbiont 16S sequence distances", "",
             "The five displayed rows are the deposited DUSA CP070274.1 reference, local DUSA-like OTU2, the close sequence relative OTU287, and the locally annotated Mesoplasma OTU3 and OTU145. Taxonomic labels follow the project; the distance calculation does not independently validate them.", "",
             f"The saved GenBank record's single annotated 16S gene is extracted directly and checked against the retained CP070274_16S.fasta. OTU2 uniquely matches positions {reference_start}–{reference_end} of the oriented 1530-bp gene. All four local sequences occupy the same 101 ACGT columns of the existing forced backbone alignment and agree with the original local OTU FASTA. The plotted DUSA row is the extracted reference fragment, with no invented bases or extension beyond the measured overlap.", "",
             "Distances are uncorrected mismatches divided by the 101 homologous ACGT sites. The compact alignment collapses only runs of at least three sites that are identical in all five sequences; the CSV and FASTA retain the complete fragment.", "",
             "| Sequence pair | Mismatches | Distance |", "|---|---:|---:|"]
    lines += [f"| {r['label_a']} / {r['label_b']} | {r['mismatches']}/101 | {r['distance_percent']:.2f}% |" for r in rows]
    lines += ["", provenance["limitation"], "",
              "Reproduce from the project root with `python3 scripts/S08_sequence_distance.py` (Pillow required). No network download is performed.", "",
              "- `figures/supplement/S08_symbiont_16S_distance.png`: final PNG, rendered directly at 300 dpi.",
              "- `output/tables/S08_symbiont_16S_distance_data.csv`: all ten unique pairwise comparisons.",
              "- `output/symbiont_16S_distance/aligned_sequences.fasta`: complete measured 101-bp fragments.",
              "- `output/symbiont_16S_distance/provenance.json`: reference coordinates, retained alignment columns, methods, limits and source SHA256 hashes.", ""]
    (DATA / "README.md").write_text("\n".join(lines))


def main():
    sequences, provenance = sequences_and_provenance()
    rows = distances(sequences, provenance)
    write_outputs(sequences, rows, provenance)
    render(sequences, rows, provenance)
    for row in rows:
        print(f"{row['label_a']} / {row['label_b']}: {row['mismatches']}/101 ({row['distance_percent']:.2f}%)")
    print(FIGURE)
    print(CSV)


if __name__ == "__main__":
    main()
