#!/usr/bin/env python3
"""Build species x family copy-number matrix for BadiRate.

Input FASTA headers should be like:
  >Species|GeneID
or
  >Species

By default, the script scans *.mnf.trimAl in the given directory.
"""

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Build count_matrix.tsv for BadiRate")
    p.add_argument("--indir", default="phylogeny_construction", help="Directory containing family FASTA files")
    p.add_argument("--pattern", default="*.mnf.trimAl", help="Glob pattern for family files")
    p.add_argument("--species-list", default=None, help="Optional species list file (one species per line)")
    p.add_argument("--out", default="count_matrix.tsv", help="Output TSV path")
    return p.parse_args()


def iter_species_ids(fasta):
    with Path(fasta).open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            hdr = line[1:].strip().split()[0]
            yield hdr.split("|", 1)[0]


def load_species(path, observed):
    if not path:
        return sorted(observed)
    species = []
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if name:
                species.append(name)
    return species


def main():
    args = parse_args()
    indir = Path(args.indir)
    families = sorted(indir.glob(args.pattern))
    if not families:
        raise SystemExit("No family files found: {}".format(indir / args.pattern))

    counts = defaultdict(lambda: defaultdict(int))
    observed_species = set()

    for fp in families:
        fam = fp.name
        if fam.endswith(".mnf.trimAl"):
            fam = fam[: -len(".mnf.trimAl")]
        for sp in iter_species_ids(fp):
            counts[sp][fam] += 1
            observed_species.add(sp)

    species = load_species(args.species_list, observed_species)
    family_names = [
        fp.name[: -len(".mnf.trimAl")] if fp.name.endswith(".mnf.trimAl") else fp.stem
        for fp in families
    ]

    out = Path(args.out)
    with out.open("w") as w:
        w.write("Species\t" + "\t".join(family_names) + "\n")
        for sp in species:
            row = [str(counts[sp].get(fam, 0)) for fam in family_names]
            w.write(sp + "\t" + "\t".join(row) + "\n")

    print("Wrote {} with {} species x {} families".format(out, len(species), len(family_names)))


if __name__ == "__main__":
    main()
