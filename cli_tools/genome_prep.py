#!/usr/bin/env python3
"""
Genome preparation for promoter/motif discovery.

Turns a genome assembly (FASTA) plus its annotation (GFF3 or GTF) into the
upstream-of-TSS sequences used by STREME, and wraps the external preprocessing
tools documented in the pipeline config so they can be run as commands.

Subcommands:
  extract-promoters  Pull N bp upstream of each gene's TSS into a FASTA
  model-repeats      One-time: build a species-specific RepeatMasker library
                     with RepeatModeler (BuildDatabase + RepeatModeler)
  mask               Repeat/low-complexity masking (RepeatMasker or dust)
  trf                Tandem Repeats Finder masking (catches (AT)n etc. that
                     TE libraries miss; commonly chained after `mask`)
  background         Markov background model (fasta-get-markov)
  run-streme         De novo motif discovery (STREME)
  run-fimo           Motif scan (FIMO) — locate STREME motif hits at controlled p/q
  run-tomtom         Motif annotation (TOMTOM) — match discovered motifs to a
                     reference TF database (PlantTFDB / JASPAR plants / CIS-BP)

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

def _require_tool(tool, executable=None):
    """
    Resolve a tool to a runnable path.

    If `executable` is given (an absolute/relative path or a command name) it is
    used instead of `tool`, so the caller can point at a module-provided binary
    that is not on PATH. Otherwise `tool` is looked up on PATH.
    """
    if executable:
        path = shutil.which(executable)
        if path is None and os.path.isfile(executable):
            path = executable  # accept an explicit file path even if not +x-flagged
        if path is None:
            raise FileNotFoundError(
                f"Specified executable for '{tool}' was not found or is not "
                f"runnable: {executable}"
            )
        return path
    path = shutil.which(tool)
    if path is None:
        raise FileNotFoundError(
            f"Required tool '{tool}' was not found on PATH. Install it (e.g. via "
            f"the MEME Suite / RepeatMasker / your HPC module), or point at it "
            f"explicitly with the corresponding --*-path option."
        )
    return path


def _run(cmd, description, cwd=None, check=True):
    print(f"\n[{description}] {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=check, cwd=cwd)


def mask_sequences(input_fasta, output=None, masker="repeatmasker",
                   species=None, threads=1, executable=None, library=None):
    """Repeat/low-complexity masking. Returns path to the masked FASTA.

    `executable` optionally overrides the path to RepeatMasker/dust.
    `library` is a custom RepeatMasker library FASTA (-lib). It takes precedence
    over `species` when both are given, since a species-specific library built
    e.g. by RepeatModeler is the recommended source of truth.
    """
    masker = masker.lower()
    if masker == "repeatmasker":
        binary = _require_tool("RepeatMasker", executable)
        if library:
            lib_path = Path(library)
            if not lib_path.exists():
                raise FileNotFoundError(
                    f"--mask-lib not found: {library} "
                    f"(resolved: {lib_path.resolve()}). "
                    f"Build it first with `model-repeats` and check the path "
                    f"(use an absolute path for safety)."
                )
            if lib_path.stat().st_size == 0:
                raise ValueError(
                    f"--mask-lib is empty: {lib_path.resolve()} — "
                    f"the RepeatModeler run probably failed; check its log."
                )
        cmd = [binary, "-pa", str(threads)]
        if library:
            cmd += ["-lib", str(library)]
            if species:
                print("[mask:RepeatMasker] both --lib and --species given; using --lib")
        elif species:
            cmd += ["-species", species]
        cmd.append(str(input_fasta))
        _run(cmd, "mask:RepeatMasker")
        produced = f"{input_fasta}.masked"
        if output and produced != str(output):
            shutil.move(produced, output)
            return str(output)
        return produced
    elif masker == "dust":
        binary = _require_tool("dust", executable)
        out = output or f"{input_fasta}.dusted.masked"
        with open(out, "w") as fh:
            subprocess.run([binary, str(input_fasta)], check=True, stdout=fh)
        print(f"[mask:dust] wrote {out}")
        return out
    else:
        raise ValueError(f"Unknown masker '{masker}' (expected 'repeatmasker' or 'dust')")


def run_trf(input_fasta, output=None, match=2, mismatch=7, delta=7,
            pm=80, pi=10, minscore=50, maxperiod=500, executable=None):
    """
    Tandem Repeats Finder masking (soft-masks tandem repeats — (AT)n, (GAA)n
    etc. — that RepeatMasker libraries typically miss). Returns the masked
    FASTA path.

    TRF deliberately exits with a non-zero status (the number of parameter
    sets it processed), so the exit code is ignored and success is verified
    by the presence of the expected `<input>.<params>.mask` output file. The
    default parameters `2 7 7 80 10 50 500` are TRF's standard recommended set.
    """
    binary = _require_tool("trf", executable)
    input_path = Path(input_fasta).resolve()
    suffix = f".{match}.{mismatch}.{delta}.{pm}.{pi}.{minscore}.{maxperiod}.mask"
    produced = input_path.parent / (input_path.name + suffix)

    cmd = [binary, str(input_path), str(match), str(mismatch), str(delta),
           str(pm), str(pi), str(minscore), str(maxperiod), "-m", "-h", "-d"]
    # TRF writes its outputs (.dat, .mask, .html) into CWD, so run from the
    # input file's directory and let it scatter sidecar files there.
    _run(cmd, "mask:TRF", cwd=str(input_path.parent), check=False)

    if not produced.exists():
        raise FileNotFoundError(
            f"TRF did not produce expected mask file: {produced}. "
            f"Verify trf ran (check stderr) and the input FASTA is non-empty."
        )
    if output and str(produced) != str(output):
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(produced, output)
        produced = Path(output)
    print(f"[mask:TRF] wrote {produced}")
    return str(produced)


def model_repeats(genome_fasta, output_dir, name=None, threads=1, ltr_struct=True,
                  legacy=False,
                  repeatmodeler_executable=None, builddatabase_executable=None):
    """
    Build a species-specific repeat library with RepeatModeler.

    Runs `BuildDatabase` then `RepeatModeler` inside `output_dir` so all their
    intermediate files (`RM_*/`, BLAST DB) stay co-located with the resulting
    library. The library FASTA is what you pass to
    `mask_sequences(..., library=...)` (RepeatMasker's -lib).

    Modes:
      - Default (RepeatModeler 2.x): writes `<name>-families.fa` and uses
        `-threads N` plus `-LTRStruct` (unless `ltr_struct=False`).
      - `legacy=True` (RepeatModeler 1.x): writes the library to
        `RM_*/consensi.fa.classified` and uses `-pa N` with `-engine ncbi`;
        `-LTRStruct` is unavailable in 1.x and is silently dropped.

    Long-running: typically hours to days for a plant genome. Designed to be
    run once per species, not on every prepare invocation.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    db_name = name or Path(genome_fasta).stem

    bd = _require_tool("BuildDatabase", builddatabase_executable)
    rm = _require_tool("RepeatModeler", repeatmodeler_executable)

    genome_abs = os.path.abspath(str(genome_fasta))
    _run([bd, "-name", db_name, genome_abs],
         "model-repeats:BuildDatabase", cwd=str(output_dir))

    if legacy:
        cmd = [rm, "-database", db_name, "-engine", "ncbi", "-pa", str(threads)]
        if ltr_struct:
            print("[model-repeats] --legacy: -LTRStruct is unavailable in "
                  "RepeatModeler 1.x; dropping it.")
    else:
        cmd = [rm, "-database", db_name, "-threads", str(threads)]
        if ltr_struct:
            cmd.append("-LTRStruct")
    _run(cmd, "model-repeats:RepeatModeler", cwd=str(output_dir))

    # 2.x writes <name>-families.fa; 1.x writes RM_*/consensi.fa.classified.
    families = output_dir / f"{db_name}-families.fa"
    if families.exists():
        library = families
    else:
        legacy_hits = sorted(output_dir.glob("RM_*/consensi.fa.classified"))
        library = legacy_hits[-1] if legacy_hits else families

    if library.exists():
        print(f"[model-repeats] custom library written -> {library}")
        if legacy and library.name != f"{db_name}-families.fa":
            print(f"[model-repeats] (pass {library} as --mask-lib)")
        return str(library)
    else:
        print(f"[model-repeats] WARNING: no library found under {output_dir} "
              f"(looked for {db_name}-families.fa and RM_*/consensi.fa.classified); "
              "check the RepeatModeler log for errors.")
        return None


def build_background(input_fasta, output="background.txt", order=1, executable=None):
    """Markov background model via fasta-get-markov (`executable` overrides path)."""
    binary = _require_tool("fasta-get-markov", executable)
    with open(output, "w") as fh:
        subprocess.run([binary, "-m", str(order), str(input_fasta)],
                       check=True, stdout=fh)
    print(f"[background] wrote {output} (order {order})")
    return output


def run_tomtom(query_motifs, target_db, output_dir, thresh=0.1,
               evalue=False, no_ssc=False, min_overlap=5, dist="pearson",
               executable=None):
    """
    TOMTOM motif-vs-database comparison: match query motifs (typically
    streme.txt) against a reference motif database (PlantTFDB, JASPAR plants,
    CIS-BP, ...) and report putative TF assignments at controlled significance.

    `thresh` defaults to 0.1 (q-value). Set `evalue=True` to interpret it as
    an E-value cutoff instead. `dist` selects the column-similarity metric:
    pearson (default), ed, kullback, sandelin, or allr.
    `executable` overrides the path to the tomtom binary.
    """
    binary = _require_tool("tomtom", executable)
    cmd = [binary, "-oc", str(output_dir), "-thresh", str(thresh),
           "-min-overlap", str(min_overlap), "-dist", dist]
    if evalue:
        cmd.append("-evalue")
    if no_ssc:
        cmd.append("-no-ssc")
    cmd += [str(query_motifs), str(target_db)]
    _run(cmd, "run-tomtom")
    return str(output_dir)


def run_fimo(motif_file, sequence_file, output_dir, thresh=1e-4,
             qv_thresh=False, no_qvalue=False, max_strand=False,
             parse_genomic_coord=False, bgfile=None, motif=None,
             executable=None):
    """
    FIMO motif scan: locate matches of `motif_file` (e.g. streme.txt, MEME format)
    in `sequence_file` (FASTA), writing results to `output_dir` (--oc).

    Standard "STREME → FIMO" workflow: scan the same promoter sequences with the
    de novo motifs to localise hits at a controlled p-/q-value. Pass `bgfile` to
    use a sequence-derived background (recommended); the literal string
    "motif-file" tells FIMO to use the background embedded in the motif file.
    `executable` overrides the path to the fimo binary.
    """
    binary = _require_tool("fimo", executable)
    cmd = [binary, "--oc", str(output_dir), "--thresh", str(thresh)]
    if qv_thresh:
        cmd.append("--qv-thresh")
    if no_qvalue:
        cmd.append("--no-qvalue")
    if max_strand:
        cmd.append("--max-strand")
    if parse_genomic_coord:
        cmd.append("--parse-genomic-coord")
    if bgfile:
        cmd += ["--bgfile", str(bgfile)]
    if motif:
        cmd += ["--motif", str(motif)]
    cmd += [str(motif_file), str(sequence_file)]
    _run(cmd, "run-fimo")
    return str(output_dir)


def run_streme(input_fasta, output_dir, nmotifs=200, minw=6, maxw=20,
               thresh=0.05, background=None, threads=1, executable=None):
    """STREME motif discovery (`executable` overrides the path to streme)."""
    binary = _require_tool("streme", executable)
    cmd = [binary, "--p", str(input_fasta),
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
    p.add_argument("--lib", dest="library",
                   help="Custom RepeatMasker library FASTA (-lib); recommended for "
                        "non-model species. Takes precedence over --species.")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--masker-path", dest="executable",
                   help="Path to the RepeatMasker/dust executable (overrides PATH lookup)")

    p = sub.add_parser("model-repeats",
                       help="Build a species-specific RepeatMasker library with RepeatModeler")
    p.add_argument("genome_fasta")
    p.add_argument("--output-dir", "-o", required=True,
                   help="Directory for the BLAST database and the resulting <name>-families.fa library")
    p.add_argument("--name", help="Database name (default: stem of the genome FASTA filename)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--no-ltr-struct", action="store_true",
                   help="Disable RepeatModeler's -LTRStruct stage (faster but misses LTR families)")
    p.add_argument("--legacy", action="store_true",
                   help="Use RepeatModeler 1.x flags: -pa instead of -threads, -engine ncbi, "
                        "no -LTRStruct. Library is written to RM_*/consensi.fa.classified.")
    p.add_argument("--repeatmodeler-path", dest="repeatmodeler_executable",
                   help="Path to the RepeatModeler executable (overrides PATH lookup)")
    p.add_argument("--builddatabase-path", dest="builddatabase_executable",
                   help="Path to the BuildDatabase executable (overrides PATH lookup)")

    p = sub.add_parser("trf", help="Tandem Repeats Finder masking (chains nicely after `mask`)")
    p.add_argument("input_fasta")
    p.add_argument("--output", "-o", help="Output masked FASTA path")
    p.add_argument("--match", type=int, default=2)
    p.add_argument("--mismatch", type=int, default=7)
    p.add_argument("--delta", type=int, default=7)
    p.add_argument("--pm", type=int, default=80, help="Match probability (default: 80)")
    p.add_argument("--pi", type=int, default=10, help="Indel probability (default: 10)")
    p.add_argument("--minscore", type=int, default=50)
    p.add_argument("--maxperiod", type=int, default=500)
    p.add_argument("--trf-path", dest="executable",
                   help="Path to the trf executable (overrides PATH lookup)")

    p = sub.add_parser("background", help="Markov background model")
    p.add_argument("input_fasta")
    p.add_argument("--output", "-o", default="background.txt")
    p.add_argument("--order", type=int, default=1)
    p.add_argument("--fasta-get-markov-path", dest="executable",
                   help="Path to the fasta-get-markov executable (overrides PATH lookup)")

    p = sub.add_parser("run-streme", help="Run STREME motif discovery")
    p.add_argument("input_fasta")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--nmotifs", type=int, default=200)
    p.add_argument("--minw", type=int, default=6)
    p.add_argument("--maxw", type=int, default=20)
    p.add_argument("--thresh", type=float, default=0.05)
    p.add_argument("--background", help="Background model file (--bfile)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--streme-path", dest="executable",
                   help="Path to the streme executable (overrides PATH lookup)")

    p = sub.add_parser("run-tomtom",
                       help="TOMTOM motif-vs-database comparison "
                            "(matches STREME motifs to reference TF PWMs)")
    p.add_argument("query_motifs", help="Query motifs (MEME format, e.g. streme.txt)")
    p.add_argument("target_db", help="Target motif database (MEME format; e.g. PlantTFDB Ath)")
    p.add_argument("--output-dir", "-o", required=True, help="TOMTOM output directory (-oc)")
    p.add_argument("--thresh", type=float, default=0.1,
                   help="Significance threshold (q-value by default; E-value with --evalue) "
                        "(default: 0.1)")
    p.add_argument("--evalue", action="store_true",
                   help="Interpret --thresh as an E-value cutoff instead of a q-value")
    p.add_argument("--no-ssc", action="store_true",
                   help="Disable small-sample correction")
    p.add_argument("--min-overlap", type=int, default=5,
                   help="Minimum overlap between query and target columns (default: 5)")
    p.add_argument("--dist", choices=["pearson", "ed", "kullback", "sandelin", "allr"],
                   default="pearson",
                   help="Column-similarity metric (default: pearson)")
    p.add_argument("--tomtom-path", dest="executable",
                   help="Path to the tomtom executable (overrides PATH lookup)")

    p = sub.add_parser("run-fimo", help="Run FIMO motif scan (typically with STREME motifs)")
    p.add_argument("motif_file", help="MEME-format motif file (e.g. streme.txt from STREME output)")
    p.add_argument("sequence_file", help="FASTA of sequences to scan")
    p.add_argument("--output-dir", "-o", required=True, help="FIMO output directory (--oc)")
    p.add_argument("--thresh", type=float, default=1e-4,
                   help="Match p-value threshold, or q-value when --qv-thresh (default: 1e-4)")
    p.add_argument("--qv-thresh", action="store_true",
                   help="Interpret --thresh as a q-value cutoff (e.g. 0.05)")
    p.add_argument("--no-qvalue", action="store_true",
                   help="Skip q-value computation (recommended when scanning very large databases)")
    p.add_argument("--max-strand", action="store_true",
                   help="Report only the higher-scoring strand of overlapping matches")
    p.add_argument("--parse-genomic-coord", action="store_true",
                   help="Parse genomic coordinates from FASTA headers (e.g. chr1:1000-2000)")
    p.add_argument("--bgfile",
                   help="Background model file. Use the literal 'motif-file' to use "
                        "the background embedded in the motif file (default).")
    p.add_argument("--motif",
                   help="Restrict to a specific motif ID from the motif file")
    p.add_argument("--fimo-path", dest="executable",
                   help="Path to the fimo executable (overrides PATH lookup)")

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
                       species=args.species, threads=args.threads,
                       executable=args.executable, library=args.library)
    elif args.command == "model-repeats":
        model_repeats(args.genome_fasta, args.output_dir,
                      name=args.name, threads=args.threads,
                      ltr_struct=not args.no_ltr_struct, legacy=args.legacy,
                      repeatmodeler_executable=args.repeatmodeler_executable,
                      builddatabase_executable=args.builddatabase_executable)
    elif args.command == "trf":
        run_trf(args.input_fasta, output=args.output,
                match=args.match, mismatch=args.mismatch, delta=args.delta,
                pm=args.pm, pi=args.pi, minscore=args.minscore,
                maxperiod=args.maxperiod, executable=args.executable)
    elif args.command == "background":
        build_background(args.input_fasta, output=args.output, order=args.order,
                         executable=args.executable)
    elif args.command == "run-streme":
        run_streme(args.input_fasta, args.output_dir, nmotifs=args.nmotifs,
                   minw=args.minw, maxw=args.maxw, thresh=args.thresh,
                   background=args.background, threads=args.threads,
                   executable=args.executable)
    elif args.command == "run-tomtom":
        run_tomtom(args.query_motifs, args.target_db, args.output_dir,
                   thresh=args.thresh, evalue=args.evalue, no_ssc=args.no_ssc,
                   min_overlap=args.min_overlap, dist=args.dist,
                   executable=args.executable)
    elif args.command == "run-fimo":
        run_fimo(args.motif_file, args.sequence_file, args.output_dir,
                 thresh=args.thresh, qv_thresh=args.qv_thresh,
                 no_qvalue=args.no_qvalue, max_strand=args.max_strand,
                 parse_genomic_coord=args.parse_genomic_coord,
                 bgfile=args.bgfile, motif=args.motif,
                 executable=args.executable)
    return 0


if __name__ == "__main__":
    sys.exit(main())
