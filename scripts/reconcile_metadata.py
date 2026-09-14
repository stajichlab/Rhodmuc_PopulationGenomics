#!/usr/bin/env python3
"""Reconcile samplesheet.csv strain IDs against metadata.txt species/ploidy calls.

Produces a strain -> confirmed_species -> ploidy source-of-truth table
(reconciliation_table.csv) plus a short text report of the specific
inconsistencies found, so the mucilaginosa-only scope for this project can be
confirmed by hand before joint genotyping.

Inputs (read-only, all relative to the project root):
    samplesheet.csv   - patient,sex,status,sample,lane,fastq_1,fastq_2
    metadata.txt      - sample,hashtag,species,strain,ploidy,MAT_type,Origin,Environment (tab-separated)
    full_samples.csv  - Strain,FileBase,Note (curated reclassification notes)
    missing_strains.tsv - Strain,MissingFiles (strains awaiting FASTQ files)

Outputs:
    reconciliation_table.csv - one row per strain in samplesheet.csv
    reconciliation_report.txt - human-readable summary of flagged issues
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

TRUE_MUCILAGINOSA = "R. mucilaginosa"
AFF_MUCILAGINOSA = "R. aff. mucilaginosa"


def load_samplesheet(path):
    strains = defaultdict(int)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            strains[row["sample"].strip()] += 1
    return strains


def load_metadata(path):
    """Return dict: bare_strain -> list of metadata rows (dicts).

    The metadata 'sample' column carries a species-prefixed name and may
    include computationally derived _hap1/_hap2 rows for phased hybrids; the
    'strain' column is the bare physical-strain id that samplesheet.csv uses.
    """
    rows_by_strain = defaultdict(list)
    all_rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            all_rows.append(row)
            rows_by_strain[row["strain"]].append(row)
    return rows_by_strain, all_rows


def load_full_samples_notes(path):
    notes = {}
    if not path.exists():
        return notes
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            notes[row["Strain"].strip()] = row.get("Note", "").strip()
    return notes


def load_missing_strains(path):
    missing = set()
    if not path.exists():
        return missing
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            missing.add(row["Strain"].strip())
    return missing


def find_duplicate_metadata_rows(all_rows):
    seen = {}
    dups = []
    for row in all_rows:
        key = tuple(row.get(c, "") for c in
                     ("sample", "hashtag", "species", "strain", "ploidy", "MAT_type", "Origin", "Environment"))
        if key in seen:
            dups.append(row["sample"])
        else:
            seen[key] = True
    return dups


def classify_species(species):
    if species == TRUE_MUCILAGINOSA:
        return "INCLUDE"
    if species == AFF_MUCILAGINOSA:
        return "BORDERLINE_AFF"
    if species == "":
        return "NO_SPECIES_CALL"
    return "EXCLUDE_OTHER_SPECIES"


def pick_primary_row(rows):
    """Prefer a diploid parent row over derived _hap1/_hap2 rows when both exist."""
    non_hap = [r for r in rows if not r["sample"].endswith(("_hap1", "_hap2"))]
    return non_hap[0] if non_hap else rows[0]


def main():
    samplesheet_strains = load_samplesheet(ROOT / "samplesheet.csv")
    meta_by_strain, all_meta_rows = load_metadata(ROOT / "metadata.txt")
    notes = load_full_samples_notes(ROOT / "full_samples.csv")
    missing = load_missing_strains(ROOT / "missing_strains.tsv")
    dup_samples = find_duplicate_metadata_rows(all_meta_rows)

    out_rows = []
    strains_missing_metadata = []
    strains_species_mismatch_name_vs_field = []

    for strain in sorted(samplesheet_strains):
        n_lanes = samplesheet_strains[strain]
        meta_rows = meta_by_strain.get(strain, [])
        note = notes.get(strain, "")
        if not meta_rows:
            strains_missing_metadata.append(strain)
            out_rows.append({
                "strain": strain,
                "n_lanes": n_lanes,
                "species": "",
                "ploidy": "",
                "MAT_type": "",
                "Origin": "",
                "Environment": "",
                "full_samples_note": note,
                "scope": "NO_METADATA",
            })
            continue

        primary = pick_primary_row(meta_rows)
        species = primary["species"]
        scope = classify_species(species)

        # Flag cases where the metadata 'sample' name still says mucilaginosa
        # but the species field says otherwise (the reclassification-lag issue).
        if "mucilaginosa" in primary["sample"].lower() and species not in (TRUE_MUCILAGINOSA, AFF_MUCILAGINOSA, ""):
            strains_species_mismatch_name_vs_field.append((strain, primary["sample"], species))

        out_rows.append({
            "strain": strain,
            "n_lanes": n_lanes,
            "species": species,
            "ploidy": primary["ploidy"],
            "MAT_type": primary["MAT_type"],
            "Origin": primary["Origin"],
            "Environment": primary["Environment"],
            "full_samples_note": note,
            "scope": scope,
        })

    # Strains that only exist in metadata.txt (no sequence data at all)
    meta_only_strains = sorted(set(meta_by_strain) - set(samplesheet_strains))

    out_csv = ROOT / "reconciliation_table.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "strain", "n_lanes", "species", "ploidy", "MAT_type",
            "Origin", "Environment", "full_samples_note", "scope",
        ])
        writer.writeheader()
        writer.writerows(out_rows)

    include_n = sum(1 for r in out_rows if r["scope"] == "INCLUDE")
    borderline_n = sum(1 for r in out_rows if r["scope"] == "BORDERLINE_AFF")
    exclude_n = sum(1 for r in out_rows if r["scope"] == "EXCLUDE_OTHER_SPECIES")
    no_meta_n = sum(1 for r in out_rows if r["scope"] == "NO_METADATA")
    no_species_n = sum(1 for r in out_rows if r["scope"] == "NO_SPECIES_CALL")

    report_lines = []
    report_lines.append("Reconciliation report: samplesheet.csv (%d strains) vs metadata.txt (%d strains)"
                         % (len(samplesheet_strains), len(meta_by_strain)))
    report_lines.append("")
    report_lines.append("Scope counts among sequenced (samplesheet.csv) strains:")
    report_lines.append(f"  INCLUDE (species == '{TRUE_MUCILAGINOSA}')        : {include_n}")
    report_lines.append(f"  BORDERLINE_AFF (species == '{AFF_MUCILAGINOSA}') : {borderline_n}")
    report_lines.append(f"  EXCLUDE_OTHER_SPECIES                            : {exclude_n}")
    report_lines.append(f"  NO_SPECIES_CALL (metadata row exists, blank species): {no_species_n}")
    report_lines.append(f"  NO_METADATA (no metadata.txt row at all)         : {no_meta_n}")
    report_lines.append("")

    report_lines.append(f"Strains in samplesheet.csv with NO metadata.txt row ({len(strains_missing_metadata)}):")
    for s in strains_missing_metadata:
        report_lines.append(f"  - {s}")
    report_lines.append("")

    report_lines.append("Strains where metadata 'sample' name says mucilaginosa but 'species' field "
                         f"disagrees ({len(strains_species_mismatch_name_vs_field)}):")
    for strain, sample_name, species in strains_species_mismatch_name_vs_field:
        report_lines.append(f"  - {strain}: sample='{sample_name}' species='{species}'")
    report_lines.append("")

    report_lines.append(f"Exact duplicate metadata.txt rows ({len(dup_samples)}):")
    for s in dup_samples:
        report_lines.append(f"  - {s}")
    report_lines.append("")

    report_lines.append(f"Strains in metadata.txt with no sequence data in samplesheet.csv ({len(meta_only_strains)}):")
    for s in meta_only_strains[:50]:
        report_lines.append(f"  - {s}")
    if len(meta_only_strains) > 50:
        report_lines.append(f"  ... and {len(meta_only_strains) - 50} more")
    report_lines.append("")

    still_missing_files = sorted(set(missing) & set(samplesheet_strains))
    report_lines.append(f"Strains flagged in missing_strains.tsv that ALSO have sequence data in samplesheet.csv "
                         f"(should be reviewed - stale entry?) ({len(still_missing_files)}):")
    for s in still_missing_files:
        report_lines.append(f"  - {s}")

    report_path = ROOT / "reconciliation_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")

    print(f"Wrote {out_csv}")
    print(f"Wrote {report_path}")
    print()
    print("\n".join(report_lines[:20]))


if __name__ == "__main__":
    main()
