#!/usr/bin/env python3
"""Cross-check genotyped strain IDs and read files against the strain DB.

For every strain in samplesheet.csv (the reads behind each CRAM):
  - resolve the strain name to a DB strain_id via strain_aliases
    (alias_normalized = upper case, non-alphanumerics removed)
  - check that every FASTQ we used is listed in sequencing_runs for that
    strain_id (filebase with R[12], or sra_run); report files that the DB
    assigns to a different strain, or does not list at all
  - compare species: DB strains.species_full, metadata.txt, the phenotype
    table, and the genetic placement (results/variant_qc/divergence/strain_placement.tsv)

Run with the StrainDB pixi python (has duckdb), from the project root:
  /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_StrainDB/.pixi/envs/default/bin/python \
      scripts/variant_qc/sample_id_crosscheck.py
Writes results/variant_qc/sample_id_crosscheck.tsv and prints a summary.
"""
import collections
import csv
import os
import re
import duckdb

DB = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_StrainDB/db/strains.duckdb"
PHENO = "/bigdata/stajichlab/shared/projects/Rhodotorula/ExRhodotorula_Phenotypes/strains.csv"
OUT = "results/variant_qc/sample_id_crosscheck.tsv"


def norm(x):
    return re.sub(r"[^A-Z0-9]", "", x.upper())


def file_key(name):
    """FASTQ basename or DB filebase -> comparable key."""
    b = os.path.basename(name)
    m = re.match(r"^(SRR\d+)(_[12])?(\.fastq\.gz)?$", b)
    if m:
        return m.group(1)
    b = b.replace("R[12]", "R1")
    b = re.sub(r"_R2_", "_R1_", b)
    return b


con = duckdb.connect(DB, read_only=True)
alias = collections.defaultdict(set)
for a, an, sid in con.execute("select alias, alias_normalized, strain_id from strain_aliases").fetchall():
    alias[an].add(sid)
    alias[norm(a)].add(sid)
for sid, in con.execute("select strain_id from strains").fetchall():
    alias[norm(sid)].add(sid)
species = dict(con.execute("select strain_id, species_full from strains").fetchall())
review = {r[0]: r[1] for r in con.execute("select strain_id, review_reason from strains where needs_review").fetchall()}
owner = collections.defaultdict(set)
db_files = collections.defaultdict(set)
for sid, sra, fb in con.execute("select strain_id, sra_run, filebase from sequencing_runs").fetchall():
    for k in filter(None, [sra] + (fb.split(";") if fb else [])):
        k = file_key(k.strip())
        owner[k].add(sid)
        db_files[sid].add(k)

ours = collections.defaultdict(list)
for r in csv.DictReader(open("samplesheet.csv")):
    ours[r["sample"]].append(file_key(r["fastq_1"]))
meta = {}
for r in csv.DictReader(open("metadata.txt"), delimiter="\t"):
    meta.setdefault(r["strain"], r["species"])
ph = {r["STRAIN"]: r["ASSIGNEDSPECIES"] for r in csv.DictReader(open(PHENO))}
place = {r["strain"]: r for r in csv.DictReader(open("results/variant_qc/divergence/strain_placement.tsv"), delimiter="\t")}


def genetic(s):
    g = place.get(s, {}).get("group", "NA")
    return {"rmuc_core": "R. mucilaginosa", "NRRL_Y-2510": "R. mucilaginosa",
            "no-metadata haploid (6)": "R. mucilaginosa"}.get(g, g)


def sp_short(x):
    return (x or "").replace("Rhodotorula", "R.").strip()


rows = []
for s in sorted(ours):
    ids = sorted(alias.get(norm(s), set()))
    sid = ids[0] if len(ids) == 1 else ""
    files = ours[s]
    st, foreign = [], []
    for k in files:
        o = owner.get(k, set())
        if not o:
            st.append("not_in_db")
        elif sid and sid in o:
            st.append("ok")
        else:
            st.append("other_strain")
            foreign.append(f"{k}->{','.join(sorted(o))}")
    file_status = "ok" if st and all(x == "ok" for x in st) else ";".join(sorted(set(st)))
    gsp = genetic(s)
    dsp = sp_short(species.get(sid, ""))
    flags = []
    if not ids:
        flags.append("name_not_in_db")
    elif len(ids) > 1:
        flags.append("name_ambiguous")
    if "other_strain" in st:
        flags.append("reads_belong_to_other_strain")
    if "not_in_db" in st:
        flags.append("reads_not_in_db")
    gen_is_rmuc = gsp == "R. mucilaginosa"
    if dsp and gen_is_rmuc and "mucilaginosa" not in dsp:
        flags.append("db_species_not_rmuc_but_genotype_rmuc")
    if dsp and gsp in ("R. frigidialcoholis", "R. aff. mucilaginosa", "candidate (no metadata)") and dsp == "R. mucilaginosa":
        flags.append("db_says_rmuc_but_genotype_other")
    p = sp_short(ph.get(s, ""))
    if p and gen_is_rmuc and "mucilaginosa" not in p:
        flags.append("phenotype_species_not_rmuc_but_genotype_rmuc")
    rows.append([s, ";".join(ids) or "NA", file_status, " ".join(foreign), str(len(files)), str(len(db_files.get(sid, ()))),
                 dsp or "NA", sp_short(meta.get(s, "")) or "NA", p or "NA", gsp, place.get(s, {}).get("group", "NA"),
                 review.get(sid, ""), ",".join(flags)])

with open(OUT, "w") as o:
    o.write("strain\tdb_strain_id\tread_files\tforeign_files\tn_files_used\tn_files_in_db\tdb_species\t"
            "metadata_species\tphenotype_species\tgenetic_species\tgenetic_group\tdb_review_reason\tflags\n")
    for r in rows:
        o.write("\t".join(r) + "\n")

print(f"strains in samplesheet: {len(rows)}")
c = collections.Counter(f for r in rows for f in r[-1].split(",") if f)
for k, v in c.most_common():
    print(f"  {k}: {v}")
print(f"  no flag: {sum(1 for r in rows if not r[-1])}")
