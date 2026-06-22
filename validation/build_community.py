#!/usr/bin/env python3
"""Build CAMISIM inputs + ground truth from genome_manifest.tsv.

Outputs (in --out):
  genome_to_id.tsv              genome_id <tab> /abs/path/to.fna   (CAMISIM)
  metadata.tsv                  genome_ID, OTU, NCBI_ID, novelty_category (CAMISIM)
  distributions/<sample>.tsv    genome_id <tab> relative_abundance (one per sample)
  truth_species.tsv             species, genus, domain, rel_abundance  (per sample)
  truth_genus.tsv               genus, domain, rel_abundance           (per sample)
  samples.txt                   list of sample names

Profiles:
  even      -> all target genomes equal
  staggered -> log-normal abundances (orders-of-magnitude spread; stresses min-reads)
Host (Homo sapiens) is added as a fixed contamination fraction and is NOT truth.
"""
import argparse, csv, math, os, random
from collections import defaultdict

def read_manifest(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for x in r:
            rows.append(x)
    return rows

def lognormal_weights(n, sigma, rng):
    w = [math.exp(rng.gauss(0.0, sigma)) for _ in range(n)]
    s = sum(w)
    return [x / s for x in w]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="genome_manifest.tsv")
    ap.add_argument("--out", default="camisim_in")
    ap.add_argument("--replicates", type=int, default=4)
    ap.add_argument("--host-frac", type=float, default=0.05,
                    help="fraction of the community that is host (contaminant)")
    ap.add_argument("--sigma", type=float, default=2.0,
                    help="log-normal sigma for the staggered profile")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    rows = read_manifest(args.manifest)
    targets = [x for x in rows if x["domain"] != "Host"]
    hosts = [x for x in rows if x["domain"] == "Host"]
    for i, x in enumerate(rows):
        x["gid"] = "G%03d" % i  # CAMISIM genome id

    os.makedirs(os.path.join(args.out, "distributions"), exist_ok=True)

    # genome_to_id + metadata
    with open(os.path.join(args.out, "genome_to_id.tsv"), "w") as f:
        for x in rows:
            f.write("%s\t%s\n" % (x["gid"], os.path.abspath(x["file"])))
    with open(os.path.join(args.out, "metadata.tsv"), "w") as f:
        f.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n")
        for i, x in enumerate(rows):
            f.write("%s\t%d\t%s\tknown_species\n" % (x["gid"], i, x["taxid"]))

    samples, truth_sp, truth_ge = [], [], []
    profiles = [("even", 0.0), ("stag", args.sigma)]
    for pname, sigma in profiles:
        for rep in range(1, args.replicates + 1):
            sample = "MOCK%s%02d" % (pname, rep)   # group = MOCKeven / MOCKstag
            samples.append(sample)
            n = len(targets)
            if sigma == 0.0:
                w = [1.0 / n] * n
            else:
                w = lognormal_weights(n, sigma, rng)
            # scale targets to (1 - host_frac), host gets host_frac
            w = [x * (1.0 - args.host_frac) for x in w]
            dist = {targets[i]["gid"]: w[i] for i in range(n)}
            for h in hosts:
                dist[h["gid"]] = args.host_frac / max(1, len(hosts))
            with open(os.path.join(args.out, "distributions", sample + ".tsv"), "w") as f:
                for gid, ab in dist.items():
                    f.write("%s\t%.10g\n" % (gid, ab))
            # truth (targets only), rolled to species + genus
            sp_ab = defaultdict(float)
            ge_ab = defaultdict(float)
            for i in range(n):
                org = targets[i]["organism"]
                dom = targets[i]["domain"]
                genus = org.split()[0]
                sp_ab[(org, genus, dom)] += w[i]
                ge_ab[(genus, dom)] += w[i]
            for (sp, ge, dom), ab in sp_ab.items():
                truth_sp.append((sample, sp, ge, dom, ab))
            for (ge, dom), ab in ge_ab.items():
                truth_ge.append((sample, ge, dom, ab))

    with open(os.path.join(args.out, "truth_species.tsv"), "w") as f:
        f.write("sample\tspecies\tgenus\tdomain\trel_abundance\n")
        for r in truth_sp:
            f.write("%s\t%s\t%s\t%s\t%.10g\n" % r)
    with open(os.path.join(args.out, "truth_genus.tsv"), "w") as f:
        f.write("sample\tgenus\tdomain\trel_abundance\n")
        for r in truth_ge:
            f.write("%s\t%s\t%s\t%.10g\n" % r)
    with open(os.path.join(args.out, "samples.txt"), "w") as f:
        f.write("\n".join(samples) + "\n")

    print("Wrote CAMISIM inputs + truth for %d samples to %s/" % (len(samples), args.out))
    print("Targets: %d genomes | Host: %d" % (len(targets), len(hosts)))

if __name__ == "__main__":
    main()
