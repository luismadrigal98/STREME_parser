#!/usr/bin/env python3
"""
Genome preparation for promoter/motif discovery.

Turns a genome assembly (FASTA) plus its annotation (GFF3 or GTF) into the
upstream-of-TSS sequences used by STREME, and wraps the external preprocessing
tools documented in the pipeline config so they can be run as commands.

Subcommands:
  extract-promoters  Pull N bp upstream of each gene's TSS into a FASTA
  mask               Repeat/low-complexity masking (RepeatMasker or dust)
  background         Markov background model (fasta-get-markov)
  run-streme         Motif discovery (STREME)

The promoter extractor is pure standard library (with an optional pyfaidx
speed-up) so it works on any machine without external bioinformatics tools.
The mask/background/run-streme steps shell out to their respective tools and
fail with a clear message when the tool is not installed.
"""

import os
import re
import sys
import shutil
import argparse
import subprocess
from pathlib import Path


# --------------------------------------------------------------------------- #
# FASTA random access (faidx-style, no full-genome load)
# --------------------------------------------------------------------------- #

class FastaIndex:
    """
    Minimal faidx-style index for random access into a (possibly large) FASTA.

    Records, per contig, the byte offset of its first base and the line layout
    so a region can be fetched by seeking instead of loading the whole genome.
    Assumes a uniform line width within each contig, which holds for standard
    genome FASTAs. Falls back to pyfaidx automatically when it is installed.
    """

    def __init__(self, fasta_path):
        self.path = str(fasta_path)
        self._pyfaidx = None
        self.records = {}  # name -> (offset, length, bases_per_line, bytes_per_line)

        try:
            import pyfaidx  # type: ignore
            self._pyfaidx = pyfaidx.Fasta(self.path, rebuild=False)
            self.records = {name: (0, len(self._pyfaidx[name]), 0, 0)
                            for name in self._pyfaidx.keys()}
            return
        except Exception:
            self._pyfaidx = None

        self._build_index()

    def _build_index(self):
        with open(self.path, "rb") as fh:
            name = None
            offset = 0
            length = 0
            bases_per_line = 0
            bytes_per_line = 0
            pos = 0
            while True:
                line = fh.readline()
                if not line:
                    break
                line_len = len(line)
                if line.startswith(b">"):
                    if name is not None:
                        self.records[name] = (offset, length, bases_per_line, bytes_per_line)
                    name = line[1:].split()[0].decode("utf-8", "replace") if len(line) > 1 else ""
                    offset = pos + line_len
                    length = 0
                    bases_per_line = 0
                    bytes_per_line = 0
                else:
                    stripped = line.rstrip(b"\r\n")
                    if bases_per_line == 0:
                        bases_per_line = len(stripped)
                        bytes_per_line = line_len
                    length += len(stripped)
                pos += line_len
            if name is not None:
                self.records[name] = (offset, length, bases_per_line, bytes_per_line)

    def __contains__(self, name):
        return name in self.records

    def contig_length(self, name):
        return self.records[name][1]

    def fetch(self, name, start, end):
        """Return uppercase sequence for [start, end) (0-based, half-open)."""
        if name not in self.records:
            raise KeyError(f"Contig '{name}' not found in {self.path}")
        offset, length, bpl, bytespl = self.records[name]
        start = max(0, start)
        end = min(length, end)
        if start >= end:
            return ""

        if self._pyfaidx is not None:
            return str(self._pyfaidx[name][start:end]).upper()

        newline_bytes = bytespl - bpl if bpl else 0
        with open(self.path, "rb") as fh:
            start_line = start // bpl
            start_col = start % bpl
            byte_start = offset + start_line * bytespl + start_col
            fh.seek(byte_start)
            # Read enough bytes to cover the region including embedded newlines.
            span = end - start
            n_newlines = (span + start_col) // bpl + 2
            to_read = span + n_newlines * max(1, newline_bytes) + bpl
            raw = fh.read(to_read)
            seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:span]
            return seq.decode("utf-8", "replace").upper()


_COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                            "TGCANtgcanYRSWMKVHDByrswmkvhdb")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def resolve_contig_filter(chromosomes):
    """
    Normalise a chromosome/contig allowlist into a set of names (or None).

    Accepts: None/empty (-> None, keep all), an iterable of names, a path to a
    file with one name per line, or a comma-separated string. Genome assemblies
    name chromosomes differently (PeChr1, Chr1, chr01, scaffolds like
    JBCEGF010000009.1), so the list is matched literally against the annotation
    seqids / FASTA names — pair it with a regex via contig_pattern when handy.
    """
    if not chromosomes:
        return None
    if isinstance(chromosomes, (set, list, tuple)):
        return {str(c).strip() for c in chromosomes if str(c).strip()}
    text = str(chromosomes)
    path = Path(text)
    if path.exists() and path.is_file():
        with open(path) as fh:
            return {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}
    return {c.strip() for c in text.split(",") if c.strip()}


# --------------------------------------------------------------------------- #
# Annotation parsing (GFF3 / GTF)
# --------------------------------------------------------------------------- #

def parse_attributes(field):
    """Parse a GFF3 ('key=value;') or GTF ('key \"value\";') attribute column."""
    attrs = {}
    field = field.strip()
    if not field:
        return attrs
    if "=" in field and '"' not in field.split("=", 1)[1][:1]:
        # GFF3 style
        for part in field.split(";"):
            part = part.strip()
            if not part or "=" not in part:
                continue
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
    else:
        # GTF style: key "value";
        for match in re.finditer(r'(\S+)\s+"([^"]*)"', field):
            attrs[match.group(1)] = match.group(2)
    return attrs


def gene_id_from_attrs(attrs, feature_type):
    """Pick the most sensible identifier from parsed attributes."""
    for key in ("ID", "gene_id", "Name", "gene_name", "locus_tag", "transcript_id"):
        if key in attrs and attrs[key]:
            value = attrs[key]
            # Strip GFF3 type prefixes like "gene:" / "mRNA:".
            if ":" in value and value.split(":", 1)[0].lower() in (
                "gene", "mrna", "transcript", "cds"
            ):
                value = value.split(":", 1)[1]
            return value
    return None


def parse_annotation(annotation_path, feature_type="gene"):
    """
    Read genes/features of interest from a GFF3 or GTF file.

    Returns a list of dicts: {gene_id, contig, start, end, strand}
    with 1-based inclusive coordinates (as in GFF/GTF).
    """
    feature_type = feature_type.lower()
    # GTF uses 'transcript'; GFF3 uses 'mRNA'. Treat them interchangeably.
    if feature_type in ("mrna", "transcript"):
        accepted = {"mrna", "transcript"}
    else:
        accepted = {feature_type}

    features = []
    seen_ids = set()
    with open(annotation_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2].lower()
            if ftype not in accepted:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            strand = cols[6] if cols[6] in ("+", "-") else "+"
            attrs = parse_attributes(cols[8])
            gid = gene_id_from_attrs(attrs, ftype)
            if gid is None:
                gid = f"{cols[0]}:{start}-{end}"
            # For multi-transcript genes keep the first occurrence per id.
            if gid in seen_ids:
                continue
            seen_ids.add(gid)
            features.append({
                "gene_id": gid,
                "contig": cols[0],
                "start": min(start, end),
                "end": max(start, end),
                "strand": strand,
            })
    return features


# --------------------------------------------------------------------------- #
# Promoter extraction
# --------------------------------------------------------------------------- #

def compute_promoter_region(feature, upstream, downstream):
    """
    Compute the genomic promoter window (1-based inclusive) for a feature.

    For '+' strand the TSS is the feature start; for '-' strand it is the
    feature end. Returns (region_start, region_end) before contig clipping.
    """
    if feature["strand"] == "-":
        tss = feature["end"]
        region_start = tss + 1 - downstream
        region_end = tss + upstream
    else:
        tss = feature["start"]
        region_start = tss - upstream
        region_end = tss - 1 + downstream
    return region_start, region_end


def _build_neighbor_bounds(features):
    """
    For --avoid-overlap: per contig, map each gene to the genomic coordinate of
    its nearest neighbour on the upstream side so the promoter can be clipped
    short of an adjacent gene's body.
    """
    by_contig = {}
    for feat in features:
        by_contig.setdefault(feat["contig"], []).append(feat)
    # left_bound[gene_id]  -> highest neighbour end strictly left of this gene
    # right_bound[gene_id] -> lowest neighbour start strictly right of this gene
    left_bound = {}
    right_bound = {}
    for contig_feats in by_contig.values():
        ordered = sorted(contig_feats, key=lambda f: (f["start"], f["end"]))
        for i, feat in enumerate(ordered):
            lb = None
            for j in range(i - 1, -1, -1):
                if ordered[j]["end"] < feat["start"]:
                    lb = ordered[j]["end"]
                    break
            rb = None
            for j in range(i + 1, len(ordered)):
                if ordered[j]["start"] > feat["end"]:
                    rb = ordered[j]["start"]
                    break
            left_bound[id(feat)] = lb
            right_bound[id(feat)] = rb
    return left_bound, right_bound


def extract_promoters(genome_fasta, annotation, output_fasta,
                      upstream=1000, downstream=0, feature_type="gene",
                      avoid_overlap=False, min_length=1, genome_name=None,
                      chromosomes=None, contig_pattern=None,
                      line_width=70, verbose=True):
    """
    Extract upstream-of-TSS promoter sequences to a FASTA file.

    Restrict to particular sequences with `chromosomes` (an allowlist; see
    resolve_contig_filter) and/or `contig_pattern` (a regex matched against the
    sequence name). A feature is kept only if it passes both filters. Returns a
    summary dict with counts.
    """
    if verbose:
        print(f"[extract-promoters] genome={genome_fasta}")
        print(f"[extract-promoters] annotation={annotation} (feature={feature_type})")
        print(f"[extract-promoters] window=-{upstream}/+{downstream} bp around TSS")

    index = FastaIndex(genome_fasta)
    features = parse_annotation(annotation, feature_type)
    if verbose:
        print(f"[extract-promoters] {len(features)} '{feature_type}' features parsed")

    # Restrict to selected chromosomes / contigs.
    allowed = resolve_contig_filter(chromosomes)
    pattern = re.compile(contig_pattern) if contig_pattern else None
    if allowed is not None or pattern is not None:
        def _keep(contig):
            if allowed is not None and contig not in allowed:
                return False
            if pattern is not None and not pattern.search(contig):
                return False
            return True
        all_contigs = sorted({f["contig"] for f in features})
        before = len(features)
        features = [f for f in features if _keep(f["contig"])]
        if verbose:
            kept = sorted({f["contig"] for f in features})
            shown = ", ".join(kept[:20]) + (" ..." if len(kept) > 20 else "")
            print(f"[extract-promoters] contig filter kept {len(features)}/{before} "
                  f"features on {len(kept)} sequence(s): {shown}")
        if not features:
            sample = ", ".join(all_contigs[:20]) + (" ..." if len(all_contigs) > 20 else "")
            print("[extract-promoters] WARNING: contig filter removed ALL features — "
                  "check --chromosomes/--contig-pattern against the annotation seqids "
                  f"(present: {sample}) and the FASTA sequence names.")

    left_bound = right_bound = None
    if avoid_overlap:
        left_bound, right_bound = _build_neighbor_bounds(features)

    written = 0
    skipped_missing_contig = 0
    skipped_short = 0
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(output_fasta, "w") as out:
        for feat in features:
            contig = feat["contig"]
            if contig not in index:
                skipped_missing_contig += 1
                continue
            contig_len = index.contig_length(contig)
            region_start, region_end = compute_promoter_region(feat, upstream, downstream)

            if avoid_overlap:
                if feat["strand"] == "-":
                    rb = right_bound.get(id(feat))
                    if rb is not None:
                        region_end = min(region_end, rb - 1)
                else:
                    lb = left_bound.get(id(feat))
                    if lb is not None:
                        region_start = max(region_start, lb + 1)

            region_start = max(1, region_start)
            region_end = min(contig_len, region_end)
            if region_end < region_start:
                skipped_short += 1
                continue

            seq = index.fetch(contig, region_start - 1, region_end)
            if feat["strand"] == "-":
                seq = reverse_complement(seq)
            if len(seq) < min_length:
                skipped_short += 1
                continue

            header = f">{feat['gene_id']}"
            if genome_name:
                header += f" genome={genome_name}"
            header += f" loc={contig}:{region_start}-{region_end}({feat['strand']})"
            out.write(header + "\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i:i + line_width] + "\n")
            written += 1

    if verbose:
        print(f"[extract-promoters] wrote {written} sequences -> {output_fasta}")
        if skipped_missing_contig:
            print(f"[extract-promoters] skipped {skipped_missing_contig} "
                  f"(contig absent from genome FASTA)")
        if skipped_short:
            print(f"[extract-promoters] skipped {skipped_short} (region < min_length)")

    return {
        "written": written,
        "skipped_missing_contig": skipped_missing_contig,
        "skipped_short": skipped_short,
        "output": output_fasta,
    }


# --------------------------------------------------------------------------- #
# External-tool wrappers (mask / background / run-streme)
# --------------------------------------------------------------------------- #

def _require_tool(tool):
    path = shutil.which(tool)
    if path is None:
        raise FileNotFoundError(
            f"Required tool '{tool}' was not found on PATH. Install it (e.g. via "
            f"the MEME Suite / RepeatMasker / your HPC module) and try again."
        )
    return path


def _run(cmd, description):
    print(f"\n[{description}] {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True)


def mask_sequences(input_fasta, output=None, masker="repeatmasker",
                   species=None, threads=1):
    """Repeat/low-complexity masking. Returns path to the masked FASTA."""
    masker = masker.lower()
    if masker == "repeatmasker":
        _require_tool("RepeatMasker")
        cmd = ["RepeatMasker", "-pa", str(threads)]
        if species:
            cmd += ["-species", species]
        cmd.append(str(input_fasta))
        _run(cmd, "mask:RepeatMasker")
        produced = f"{input_fasta}.masked"
        if output and produced != str(output):
            shutil.move(produced, output)
            return str(output)
        return produced
    elif masker == "dust":
        _require_tool("dust")
        out = output or f"{input_fasta}.dusted.masked"
        with open(out, "w") as fh:
            subprocess.run(["dust", str(input_fasta)], check=True, stdout=fh)
        print(f"[mask:dust] wrote {out}")
        return out
    else:
        raise ValueError(f"Unknown masker '{masker}' (expected 'repeatmasker' or 'dust')")


def build_background(input_fasta, output="background.txt", order=1):
    """Markov background model via fasta-get-markov."""
    _require_tool("fasta-get-markov")
    with open(output, "w") as fh:
        subprocess.run(["fasta-get-markov", "-m", str(order), str(input_fasta)],
                       check=True, stdout=fh)
    print(f"[background] wrote {output} (order {order})")
    return output


def run_streme(input_fasta, output_dir, nmotifs=200, minw=6, maxw=20,
               thresh=0.05, background=None, threads=1):
    """STREME motif discovery."""
    _require_tool("streme")
    cmd = ["streme", "--p", str(input_fasta),
           "--nmotifs", str(nmotifs),
           "--minw", str(minw), "--maxw", str(maxw),
           "--thresh", str(thresh),
           "-o", str(output_dir)]
    if background:
        cmd += ["--bfile", str(background)]
    if threads and threads > 1:
        # STREME parallelises internally; expose the knob via OMP threads.
        os.environ.setdefault("OMP_NUM_THREADS", str(threads))
    _run(cmd, "run-streme")
    return str(output_dir)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_parser():
    parser = argparse.ArgumentParser(
        description="Genome preparation: promoter extraction, masking, "
                    "background model, and STREME motif discovery.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("extract-promoters",
                       help="Extract N bp upstream of each gene's TSS")
    p.add_argument("genome_fasta", help="Genome assembly FASTA")
    p.add_argument("annotation", help="Annotation file (GFF3 or GTF)")
    p.add_argument("--output", "-o", required=True, help="Output FASTA path")
    p.add_argument("--upstream", "-u", type=int, default=1000,
                   help="Bases upstream of TSS (default: 1000)")
    p.add_argument("--downstream", "-d", type=int, default=0,
                   help="Bases downstream of TSS to include (default: 0)")
    p.add_argument("--feature-type", default="gene",
                   help="Annotation feature type to use (default: gene; "
                        "also mRNA/transcript)")
    p.add_argument("--avoid-overlap", action="store_true",
                   help="Clip the window short of an adjacent gene's body")
    p.add_argument("--min-length", type=int, default=1,
                   help="Drop promoters shorter than this (default: 1)")
    p.add_argument("--genome-name", help="Genome label to embed in FASTA headers")
    p.add_argument("--chromosomes",
                   help="Restrict to these sequences: comma-separated names "
                        "(e.g. PeChr1,PeChr2,...) or a file with one name per line")
    p.add_argument("--contig-pattern",
                   help=r"Restrict to sequences whose name matches this regex "
                        r"(e.g. '^PeChr' to keep chromosomes, drop scaffolds)")

    p = sub.add_parser("mask", help="Repeat/low-complexity masking")
    p.add_argument("input_fasta")
    p.add_argument("--output", "-o", help="Output masked FASTA path")
    p.add_argument("--masker", choices=["repeatmasker", "dust"], default="repeatmasker")
    p.add_argument("--species", help="Species for RepeatMasker")
    p.add_argument("--threads", type=int, default=1)

    p = sub.add_parser("background", help="Markov background model")
    p.add_argument("input_fasta")
    p.add_argument("--output", "-o", default="background.txt")
    p.add_argument("--order", type=int, default=1)

    p = sub.add_parser("run-streme", help="Run STREME motif discovery")
    p.add_argument("input_fasta")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--nmotifs", type=int, default=200)
    p.add_argument("--minw", type=int, default=6)
    p.add_argument("--maxw", type=int, default=20)
    p.add_argument("--thresh", type=float, default=0.05)
    p.add_argument("--background", help="Background model file (--bfile)")
    p.add_argument("--threads", type=int, default=1)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 1

    if args.command == "extract-promoters":
        extract_promoters(
            args.genome_fasta, args.annotation, args.output,
            upstream=args.upstream, downstream=args.downstream,
            feature_type=args.feature_type, avoid_overlap=args.avoid_overlap,
            min_length=args.min_length, genome_name=args.genome_name,
            chromosomes=args.chromosomes, contig_pattern=args.contig_pattern,
        )
    elif args.command == "mask":
        mask_sequences(args.input_fasta, output=args.output, masker=args.masker,
                       species=args.species, threads=args.threads)
    elif args.command == "background":
        build_background(args.input_fasta, output=args.output, order=args.order)
    elif args.command == "run-streme":
        run_streme(args.input_fasta, args.output_dir, nmotifs=args.nmotifs,
                   minw=args.minw, maxw=args.maxw, thresh=args.thresh,
                   background=args.background, threads=args.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
