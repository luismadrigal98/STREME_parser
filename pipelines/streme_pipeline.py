#!/usr/bin/env python3
"""
STREME Analysis Pipeline - Main Pipeline Steps:
  0. prepare          - Genome(s) -> promoters -> mask -> background -> STREME
  1. consolidate      - Consolidate STREME motifs across genomes
  2. validate         - Validate motif consolidation quality
  3. analyze          - Run motif-expression analysis (absolute or relative)
  4. full             - Run complete pipeline (prepare optional, then 1-3)

This orchestrates the complete pipeline from genome preparation and STREME
motif discovery through consolidation, validation, and expression analysis.
"""

import os
import csv
import sys
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Make the genome-prep helpers importable regardless of working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "cli_tools"))


def read_genome_manifest(manifest_path):
    """
    Parse a tab-separated manifest of genome/annotation pairs.

    Required columns: genome, fasta, annotation
    Optional columns: expression, chromosomes, contig_pattern, lib, species
    (`chromosomes`/`contig_pattern` are per-genome sequence filters that
    override the run-wide --chromosomes/--contig-pattern; `lib` and `species`
    give per-genome RepeatMasker library / species choices that override
    --mask-lib / --species. `lib` takes precedence over `species` when both
    end up set, since a species-specific library — e.g. one built by
    `model-repeats` (RepeatModeler) — is preferred for non-model organisms.)
    """
    def opt(row, cols, name):
        return (row[cols[name]].strip()
                if name in cols and row.get(cols[name]) else None)

    specs = []
    with open(manifest_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        for required in ("genome", "fasta", "annotation"):
            if required not in cols:
                raise ValueError(
                    f"Manifest {manifest_path} is missing required column '{required}'. "
                    f"Found: {reader.fieldnames}"
                )
        for row in reader:
            if not row.get(cols["genome"]):
                continue
            specs.append({
                "genome": row[cols["genome"]].strip(),
                "fasta": row[cols["fasta"]].strip(),
                "annotation": row[cols["annotation"]].strip(),
                "expression": opt(row, cols, "expression"),
                "chromosomes": opt(row, cols, "chromosomes"),
                "contig_pattern": opt(row, cols, "contig_pattern"),
                "lib": opt(row, cols, "lib"),
                "species": opt(row, cols, "species"),
            })
    if not specs:
        raise ValueError(f"No genome rows found in manifest {manifest_path}")
    return specs


def prepare_one_genome(spec, opts):
    """
    Run the preparation chain for a single genome:
      extract promoters -> (mask) -> (background) -> (STREME).

    Returns a result dict. Designed to be picklable / runnable in a worker
    process. Per-genome failures are captured rather than raised so a batch
    can continue.
    """
    # Re-establish the import path inside spawned workers.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "cli_tools"))
    import genome_prep

    genome = spec["genome"]
    out_root = Path(opts["output"])
    work_dir = out_root / f"{genome}_prep"
    work_dir.mkdir(parents=True, exist_ok=True)

    result = {"genome": genome, "status": "ok", "steps": {}, "error": None}
    try:
        # Per-genome chromosome filters override the run-wide defaults as a unit:
        # if a manifest row sets either field, the run-wide filters are ignored
        # for this genome (so the two filter types are never mixed across scopes).
        if spec.get("chromosomes") or spec.get("contig_pattern"):
            chromosomes = spec.get("chromosomes")
            contig_pattern = spec.get("contig_pattern")
        else:
            chromosomes = opts["chromosomes"]
            contig_pattern = opts["contig_pattern"]

        promoters = work_dir / f"{genome}_promoters.fasta"
        genome_prep.extract_promoters(
            spec["fasta"], spec["annotation"], str(promoters),
            upstream=opts["upstream"], downstream=opts["downstream"],
            feature_type=opts["feature_type"], avoid_overlap=opts["avoid_overlap"],
            min_length=opts["min_length"], genome_name=genome,
            chromosomes=chromosomes, contig_pattern=contig_pattern, verbose=True,
        )
        result["steps"]["promoters"] = str(promoters)
        streme_input = promoters

        if opts["mask"] != "none":
            # Per-genome RepeatMasker library / species override the run-wide flags
            # (the wrapper prefers library over species when both are set).
            genome_lib = spec.get("lib") or opts["mask_lib"]
            genome_species = spec.get("species") or opts["species"]
            masked = genome_prep.mask_sequences(
                str(streme_input),
                output=str(work_dir / f"{genome}_promoters.masked"),
                masker=opts["mask"], species=genome_species, threads=opts["threads"],
                executable=opts["masker_path"], library=genome_lib,
            )
            result["steps"]["masked"] = masked
            streme_input = Path(masked)

        # Optional chained masking for simple-sequence / tandem repeats that
        # primary repeat libraries miss (the typical (AT)n / (GAA)n problem).
        extra = (opts["extra_mask"] or "none").lower()
        if extra in ("dust", "both"):
            dusted = genome_prep.mask_sequences(
                str(streme_input),
                output=str(work_dir / f"{genome}_promoters.masked.dust"),
                masker="dust", threads=opts["threads"],
                executable=opts["masker_path"],
            )
            result["steps"]["extra_mask_dust"] = dusted
            streme_input = Path(dusted)
        if extra in ("trf", "both"):
            trfed = genome_prep.run_trf(
                str(streme_input),
                output=str(work_dir / f"{genome}_promoters.masked.trf"),
                executable=opts["trf_path"],
            )
            result["steps"]["extra_mask_trf"] = trfed
            streme_input = Path(trfed)

        background = None
        if opts["background"]:
            background = genome_prep.build_background(
                str(streme_input), output=str(work_dir / "background.txt"),
                order=opts["background_order"], executable=opts["markov_path"],
            )
            result["steps"]["background"] = background

        if opts["run_streme"]:
            streme_dir = out_root / f"streme_{genome}"
            genome_prep.run_streme(
                str(streme_input), str(streme_dir),
                nmotifs=opts["nmotifs"], minw=opts["minw"], maxw=opts["maxw"],
                thresh=opts["thresh"], background=background, threads=opts["threads"],
                executable=opts["streme_path"],
            )
            result["steps"]["streme"] = str(streme_dir)

            if opts["run_fimo"]:
                motif_file = streme_dir / "streme.txt"
                if not motif_file.exists():
                    raise FileNotFoundError(
                        f"FIMO requested but STREME motif file not found: {motif_file}"
                    )
                fimo_bg = background if background else "motif-file"
                fimo_dir = out_root / f"fimo_{genome}"
                genome_prep.run_fimo(
                    str(motif_file), str(streme_input), str(fimo_dir),
                    thresh=opts["fimo_thresh"],
                    qv_thresh=opts["fimo_qv_thresh"],
                    no_qvalue=opts["fimo_no_qvalue"],
                    max_strand=opts["fimo_max_strand"],
                    bgfile=fimo_bg, motif=opts["fimo_motif"],
                    executable=opts["fimo_path"],
                )
                result["steps"]["fimo"] = str(fimo_dir)

            if opts["run_tomtom"]:
                motif_file = streme_dir / "streme.txt"
                if not motif_file.exists():
                    raise FileNotFoundError(
                        f"TOMTOM requested but STREME motif file not found: {motif_file}"
                    )
                if not opts["tomtom_db"]:
                    raise ValueError("--run-tomtom requires --tomtom-db <motif database>")
                tomtom_dir = out_root / f"tomtom_{genome}"
                genome_prep.run_tomtom(
                    str(motif_file), opts["tomtom_db"], str(tomtom_dir),
                    thresh=opts["tomtom_thresh"], evalue=opts["tomtom_evalue"],
                    no_ssc=opts["tomtom_no_ssc"], min_overlap=opts["tomtom_min_overlap"],
                    dist=opts["tomtom_dist"], executable=opts["tomtom_path"],
                )
                result["steps"]["tomtom"] = str(tomtom_dir)
    except Exception as exc:  # noqa: BLE001 - surface per-genome failure to caller
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def run_prepare(args):
    """Drive genome preparation for one or many genomes, parallel across genomes."""
    if args.run_fimo and args.no_streme:
        print("Error: --run-fimo requires STREME output; remove --no-streme or run the "
              "'scan' subcommand against an existing STREME directory.")
        return False

    if args.run_tomtom:
        if args.no_streme:
            print("Error: --run-tomtom requires STREME output; remove --no-streme or run "
                  "the 'annotate' subcommand against an existing STREME directory.")
            return False
        if not args.tomtom_db:
            print("Error: --run-tomtom requires --tomtom-db <reference motif database>")
            return False
        if not Path(args.tomtom_db).is_file():
            print(f"Error: --tomtom-db file not found: {args.tomtom_db}")
            return False

    if args.manifest:
        specs = read_genome_manifest(args.manifest)
    else:
        if not (args.genome and args.fasta and args.annotation):
            print("Error: provide --manifest OR all of --genome, --fasta, --annotation")
            return False
        specs = [{
            "genome": args.genome, "fasta": args.fasta,
            "annotation": args.annotation, "expression": args.expression,
            "chromosomes": None, "contig_pattern": None,
            "lib": None, "species": None,
        }]

    opts = {
        "output": args.output,
        "upstream": args.upstream, "downstream": args.downstream,
        "feature_type": args.feature_type, "avoid_overlap": args.avoid_overlap,
        "min_length": args.min_length,
        "chromosomes": args.chromosomes, "contig_pattern": args.contig_pattern,
        "mask": args.mask, "species": args.species, "mask_lib": args.mask_lib,
        "extra_mask": args.extra_mask, "trf_path": args.trf_path,
        "masker_path": args.masker_path, "markov_path": args.fasta_get_markov_path,
        "streme_path": args.streme_path,
        "run_fimo": args.run_fimo, "fimo_thresh": args.fimo_thresh,
        "fimo_qv_thresh": args.fimo_qv_thresh, "fimo_no_qvalue": args.fimo_no_qvalue,
        "fimo_max_strand": args.fimo_max_strand, "fimo_motif": args.fimo_motif,
        "fimo_path": args.fimo_path,
        "run_tomtom": args.run_tomtom, "tomtom_db": args.tomtom_db,
        "tomtom_thresh": args.tomtom_thresh, "tomtom_evalue": args.tomtom_evalue,
        "tomtom_no_ssc": args.tomtom_no_ssc,
        "tomtom_min_overlap": args.tomtom_min_overlap,
        "tomtom_dist": args.tomtom_dist, "tomtom_path": args.tomtom_path,
        "background": not args.no_background, "background_order": args.background_order,
        "run_streme": not args.no_streme,
        "nmotifs": args.nmotifs, "minw": args.minw, "maxw": args.maxw, "thresh": args.thresh,
        "threads": args.threads,
    }
    Path(args.output).mkdir(parents=True, exist_ok=True)

    print(f"\n=== PREPARE: {len(specs)} genome(s) ===")
    print(f"Output root: {args.output}")
    print(f"Parallel genomes: {args.jobs} | threads per heavy step: {args.threads}")
    for s in specs:
        print(f"  - {s['genome']}: {s['fasta']} + {s['annotation']}")

    results = []
    if args.jobs > 1 and len(specs) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(prepare_one_genome, s, opts): s["genome"] for s in specs}
            for fut in as_completed(futures):
                results.append(fut.result())
    else:
        for s in specs:
            results.append(prepare_one_genome(s, opts))

    print("\n=== PREPARE SUMMARY ===")
    ok = 0
    for r in sorted(results, key=lambda x: x["genome"]):
        if r["status"] == "ok":
            ok += 1
            done = ", ".join(r["steps"].keys())
            print(f"  ✅ {r['genome']}: {done}")
        else:
            print(f"  ❌ {r['genome']}: {r['error']}")
    print(f"\n{ok}/{len(results)} genome(s) prepared successfully.")
    if ok and not args.no_streme:
        print(f"Next: python {sys.argv[0]} consolidate {args.output} --output outputs/consolidated_streme_sites")
    return ok == len(results)


def scan_one_genome(spec, opts):
    """
    Run FIMO for one genome (used as a worker by run_scan).

    spec: {"genome", "streme_dir", "sequence_file", "bgfile"}
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "cli_tools"))
    import genome_prep

    genome = spec["genome"]
    out_root = Path(opts["output"])
    out_root.mkdir(parents=True, exist_ok=True)
    fimo_dir = out_root / f"fimo_{genome}"

    result = {"genome": genome, "status": "ok", "output": str(fimo_dir), "error": None}
    try:
        genome_prep.run_fimo(
            spec["motif_file"], spec["sequence_file"], str(fimo_dir),
            thresh=opts["thresh"], qv_thresh=opts["qv_thresh"],
            no_qvalue=opts["no_qvalue"], max_strand=opts["max_strand"],
            parse_genomic_coord=opts["parse_genomic_coord"],
            bgfile=spec["bgfile"], motif=opts["motif"],
            executable=opts["fimo_path"],
        )
    except Exception as exc:  # noqa: BLE001
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def _resolve_scan_specs(prepared_dir, sequence_override, bgfile_override, target_genomes):
    """
    Walk a prepared/ directory and resolve per-genome FIMO inputs.

    For each `streme_<genome>/` subdir: locate streme.txt, choose a sequence
    file (--sequence override > prepared/<genome>_prep/<genome>_promoters.masked
    > the .fasta), and pick a background (--bgfile override > prepared/
    <genome>_prep/background.txt > literal 'motif-file').
    """
    prepared = Path(prepared_dir)
    if not prepared.is_dir():
        raise FileNotFoundError(f"Prepared directory not found: {prepared}")

    specs = []
    missing = []
    for entry in sorted(prepared.iterdir()):
        if not (entry.is_dir() and entry.name.startswith("streme_")):
            continue
        genome = entry.name[len("streme_"):]
        if target_genomes and genome not in target_genomes:
            continue

        motif_file = entry / "streme.txt"
        if not motif_file.exists():
            missing.append(f"{entry.name}: no streme.txt")
            continue

        if sequence_override:
            sequence_file = Path(sequence_override)
        else:
            work_dir = prepared / f"{genome}_prep"
            candidates = [
                work_dir / f"{genome}_promoters.masked",
                work_dir / f"{genome}_promoters.fasta",
            ]
            sequence_file = next((c for c in candidates if c.exists()), None)
            if sequence_file is None:
                missing.append(f"{genome}: no sequence file (looked in {work_dir})")
                continue

        if bgfile_override:
            bgfile = bgfile_override
        else:
            bg_path = prepared / f"{genome}_prep" / "background.txt"
            bgfile = str(bg_path) if bg_path.exists() else "motif-file"

        specs.append({
            "genome": genome,
            "motif_file": str(motif_file),
            "sequence_file": str(sequence_file),
            "bgfile": bgfile,
        })
    return specs, missing


def run_scan(args):
    """FIMO-scan every prepared genome's STREME motifs, in parallel across genomes."""
    target = set(args.genomes.split(",")) if args.genomes else None
    specs, missing = _resolve_scan_specs(args.prepared_dir, args.sequence,
                                         args.bgfile, target)

    if not specs:
        print("Error: no FIMO-able genomes found in", args.prepared_dir)
        for m in missing:
            print(" -", m)
        return False

    output_root = args.output or args.prepared_dir
    opts = {
        "output": output_root,
        "thresh": args.thresh, "qv_thresh": args.qv_thresh,
        "no_qvalue": args.no_qvalue, "max_strand": args.max_strand,
        "parse_genomic_coord": args.parse_genomic_coord,
        "motif": args.motif, "fimo_path": args.fimo_path,
    }
    Path(output_root).mkdir(parents=True, exist_ok=True)

    print(f"\n=== SCAN: {len(specs)} genome(s) -> {output_root} ===")
    for s in specs:
        print(f"  - {s['genome']}: motifs={s['motif_file']} seq={s['sequence_file']} bg={s['bgfile']}")
    for m in missing:
        print(f"  (skipped) {m}")

    results = []
    if args.jobs > 1 and len(specs) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(scan_one_genome, s, opts): s["genome"] for s in specs}
            for fut in as_completed(futures):
                results.append(fut.result())
    else:
        for s in specs:
            results.append(scan_one_genome(s, opts))

    print("\n=== SCAN SUMMARY ===")
    ok = 0
    for r in sorted(results, key=lambda x: x["genome"]):
        if r["status"] == "ok":
            ok += 1
            print(f"  ✅ {r['genome']}: {r['output']}")
        else:
            print(f"  ❌ {r['genome']}: {r['error']}")
    print(f"\n{ok}/{len(results)} genome(s) scanned successfully.")
    return ok == len(results)


def annotate_one_genome(spec, opts):
    """
    Run TOMTOM for one genome (worker used by run_annotate).

    spec: {"genome", "motif_file"}
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "cli_tools"))
    import genome_prep

    genome = spec["genome"]
    out_root = Path(opts["output"])
    out_root.mkdir(parents=True, exist_ok=True)
    tomtom_dir = out_root / f"tomtom_{genome}"

    result = {"genome": genome, "status": "ok", "output": str(tomtom_dir), "error": None}
    try:
        genome_prep.run_tomtom(
            spec["motif_file"], opts["target_db"], str(tomtom_dir),
            thresh=opts["thresh"], evalue=opts["evalue"], no_ssc=opts["no_ssc"],
            min_overlap=opts["min_overlap"], dist=opts["dist"],
            executable=opts["tomtom_path"],
        )
    except Exception as exc:  # noqa: BLE001
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def run_annotate(args):
    """TOMTOM-match every prepared genome's STREME motifs against a TF database."""
    prepared = Path(args.prepared_dir)
    if not prepared.is_dir():
        print(f"Error: prepared directory not found: {prepared}")
        return False
    if not Path(args.target_db).is_file():
        print(f"Error: target motif database not found: {args.target_db}")
        return False

    target = set(args.genomes.split(",")) if args.genomes else None
    specs = []
    missing = []
    for entry in sorted(prepared.iterdir()):
        if not (entry.is_dir() and entry.name.startswith("streme_")):
            continue
        genome = entry.name[len("streme_"):]
        if target and genome not in target:
            continue
        motif_file = entry / "streme.txt"
        if not motif_file.exists():
            missing.append(f"{entry.name}: no streme.txt")
            continue
        specs.append({"genome": genome, "motif_file": str(motif_file)})

    if not specs:
        print("Error: no annotate-able genomes found in", args.prepared_dir)
        for m in missing:
            print(" -", m)
        return False

    output_root = args.output or args.prepared_dir
    opts = {
        "output": output_root, "target_db": args.target_db,
        "thresh": args.thresh, "evalue": args.evalue, "no_ssc": args.no_ssc,
        "min_overlap": args.min_overlap, "dist": args.dist,
        "tomtom_path": args.tomtom_path,
    }
    Path(output_root).mkdir(parents=True, exist_ok=True)

    print(f"\n=== ANNOTATE: {len(specs)} genome(s) -> {output_root} ===")
    print(f"Target DB: {args.target_db}")
    for s in specs:
        print(f"  - {s['genome']}: {s['motif_file']}")
    for m in missing:
        print(f"  (skipped) {m}")

    results = []
    if args.jobs > 1 and len(specs) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(annotate_one_genome, s, opts): s["genome"] for s in specs}
            for fut in as_completed(futures):
                results.append(fut.result())
    else:
        for s in specs:
            results.append(annotate_one_genome(s, opts))

    print("\n=== ANNOTATE SUMMARY ===")
    ok = 0
    for r in sorted(results, key=lambda x: x["genome"]):
        if r["status"] == "ok":
            ok += 1
            print(f"  ✅ {r['genome']}: {r['output']}")
        else:
            print(f"  ❌ {r['genome']}: {r['error']}")
    print(f"\n{ok}/{len(results)} genome(s) annotated successfully.")
    return ok == len(results)


def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n{description}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with error code {e.returncode}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='STREME Analysis Pipeline - Complete motif discovery, validation, and expression analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline Steps:
  1. consolidate     - Consolidate STREME motifs across lines
  2. validate        - Validate motif consolidation quality
  3. analyze         - Run motif-expression analysis (absolute or relative)
  4. gene-specific   - Run gene-focused analysis with all motifs
  5. full            - Run complete pipeline (steps 1-3)

Examples:
  # Step 1: Consolidate motifs
  %(prog)s consolidate /path/to/streme/results --output outputs/
  
  # Step 2: Validate consolidation quality
  %(prog)s validate outputs/consolidated_streme_sites.tsv
  
  # Step 3: Expression analysis
  %(prog)s analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --output results/
  
  # Step 3: Expression analysis with custom reference
  %(prog)s analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --reference-line IM500 --output results/
  
  # Step 4: Gene-specific analysis  
  %(prog)s gene-specific outputs/consolidated_streme_sites.tsv expression.tsv --genes AT1G01010,AT1G01020 --type relative
  
  # Full pipeline
  %(prog)s full /path/to/streme/results expression.tsv --output outputs/ --analysis-type relative

Expression Analysis Types:
  - absolute: Per-line analysis (each line analyzed independently)  
  - relative: Comparative analysis (each line compared to reference baseline)
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline steps')

    # Prepare subcommand: genome(s) -> promoters -> mask -> background -> STREME
    prepare_parser = subparsers.add_parser(
        'prepare',
        help='Prepare genome(s): extract promoters, mask, background, run STREME',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single genome (parallelise the heavy steps with --threads)
  %(prog)s prepare --genome P_virgatus --fasta P_virgatus.fa \\
      --annotation P_virgatus.gff3 --upstream 1000 --threads 10 --output prepared/

  # Many genomes from a manifest, run in parallel (--jobs genomes at a time)
  %(prog)s prepare --manifest genomes.tsv --jobs 4 --threads 8 --output prepared/

  # Just extract promoters locally (no external tools needed)
  %(prog)s prepare --manifest genomes.tsv --mask none --no-background --no-streme

Manifest format (tab-separated, header required):
  genome      fasta                 annotation             expression
  P_virgatus  /data/virgatus.fa     /data/virgatus.gff3    /data/virgatus_expr.tsv
  P_barbatus  /data/barbatus.fa     /data/barbatus.gff3
        """
    )
    prepare_parser.add_argument('--manifest', help='TSV of genome/fasta/annotation[/expression] rows')
    prepare_parser.add_argument('--genome', help='Genome name (single-genome mode)')
    prepare_parser.add_argument('--fasta', help='Genome FASTA (single-genome mode)')
    prepare_parser.add_argument('--annotation', help='GFF3/GTF annotation (single-genome mode)')
    prepare_parser.add_argument('--expression', help='Optional expression table (single-genome mode)')
    prepare_parser.add_argument('--output', '-o', default='prepared/', help='Output root directory')
    prepare_parser.add_argument('--jobs', '-j', type=int, default=1,
                                help='Number of genomes to prepare in parallel (default: 1)')
    prepare_parser.add_argument('--threads', type=int, default=1,
                                help='Threads per heavy step (RepeatMasker/STREME) (default: 1)')
    # Promoter extraction
    prepare_parser.add_argument('--upstream', '-u', type=int, default=1000,
                                help='Bases upstream of TSS (default: 1000)')
    prepare_parser.add_argument('--downstream', '-d', type=int, default=0,
                                help='Bases downstream of TSS to include (default: 0)')
    prepare_parser.add_argument('--feature-type', default='gene',
                                help='Annotation feature type (default: gene)')
    prepare_parser.add_argument('--avoid-overlap', action='store_true',
                                help='Clip promoter windows short of adjacent genes')
    prepare_parser.add_argument('--min-length', type=int, default=1,
                                help='Drop promoters shorter than this (default: 1)')
    prepare_parser.add_argument('--chromosomes',
                                help='Restrict to these sequences: comma-separated names '
                                     '(e.g. PeChr1,PeChr2,...) or a file with one name per line. '
                                     'Per-genome overrides may be set in the manifest.')
    prepare_parser.add_argument('--contig-pattern',
                                help=r"Restrict to sequences whose name matches this regex "
                                     r"(e.g. '^PeChr' to keep chromosomes, drop scaffolds)")
    # Masking / background / STREME
    prepare_parser.add_argument('--mask', choices=['none', 'repeatmasker', 'dust'],
                                default='repeatmasker', help='Masking step (default: repeatmasker)')
    prepare_parser.add_argument('--species', help='Species for RepeatMasker (Dfam library lookup)')
    prepare_parser.add_argument('--mask-lib',
                                help='Custom RepeatMasker library FASTA (-lib); recommended for '
                                     'non-model organisms. Build one once with the `model-repeats` '
                                     'subcommand. Per-genome overrides via the manifest `lib` column.')
    prepare_parser.add_argument('--extra-mask', choices=['none', 'dust', 'trf', 'both'],
                                default='none',
                                help='After the primary masker, chain dust and/or Tandem Repeats '
                                     'Finder (TRF) to catch (AT)n / SSR tracts a TE library misses. '
                                     '"both" runs dust then trf in sequence (default: none).')
    prepare_parser.add_argument('--trf-path',
                                help='Path to the trf executable (overrides PATH lookup)')
    prepare_parser.add_argument('--masker-path',
                                help='Path to the RepeatMasker/dust executable (overrides PATH lookup)')
    prepare_parser.add_argument('--fasta-get-markov-path',
                                help='Path to the fasta-get-markov executable (overrides PATH lookup)')
    prepare_parser.add_argument('--streme-path',
                                help='Path to the streme executable (overrides PATH lookup)')
    # Optional inline FIMO scan with the STREME motifs
    prepare_parser.add_argument('--run-fimo', action='store_true',
                                help='After STREME, run FIMO to scan the same masked promoters '
                                     'with the discovered motifs (writes fimo_<genome>/)')
    prepare_parser.add_argument('--fimo-thresh', type=float, default=1e-4,
                                help='FIMO match p-value threshold, or q-value with --fimo-qv-thresh '
                                     '(default: 1e-4)')
    prepare_parser.add_argument('--fimo-qv-thresh', action='store_true',
                                help='Interpret --fimo-thresh as a q-value cutoff (e.g. 0.05)')
    prepare_parser.add_argument('--fimo-no-qvalue', action='store_true',
                                help='Skip FIMO q-value computation')
    prepare_parser.add_argument('--fimo-max-strand', action='store_true',
                                help='FIMO: report only the higher-scoring strand of overlapping matches')
    prepare_parser.add_argument('--fimo-motif',
                                help='FIMO: restrict scan to a specific motif ID')
    prepare_parser.add_argument('--fimo-path',
                                help='Path to the fimo executable (overrides PATH lookup)')
    # Optional inline TOMTOM annotation against a reference TF motif database
    prepare_parser.add_argument('--run-tomtom', action='store_true',
                                help='After STREME, run TOMTOM to match discovered motifs '
                                     'to a reference TF database (writes tomtom_<genome>/)')
    prepare_parser.add_argument('--tomtom-db',
                                help='Reference motif database for TOMTOM (MEME format), '
                                     'e.g. PlantTFDB Arabidopsis, JASPAR plants, CIS-BP')
    prepare_parser.add_argument('--tomtom-thresh', type=float, default=0.1,
                                help='TOMTOM significance threshold (q-value by default; '
                                     'E-value with --tomtom-evalue) (default: 0.1)')
    prepare_parser.add_argument('--tomtom-evalue', action='store_true',
                                help='Interpret --tomtom-thresh as an E-value cutoff')
    prepare_parser.add_argument('--tomtom-no-ssc', action='store_true',
                                help='Disable TOMTOM small-sample correction')
    prepare_parser.add_argument('--tomtom-min-overlap', type=int, default=5,
                                help='TOMTOM minimum column overlap (default: 5)')
    prepare_parser.add_argument('--tomtom-dist',
                                choices=['pearson', 'ed', 'kullback', 'sandelin', 'allr'],
                                default='pearson',
                                help='TOMTOM column-similarity metric (default: pearson)')
    prepare_parser.add_argument('--tomtom-path',
                                help='Path to the tomtom executable (overrides PATH lookup)')
    prepare_parser.add_argument('--no-background', action='store_true',
                                help='Skip the Markov background model step')
    prepare_parser.add_argument('--background-order', type=int, default=1,
                                help='Markov order for background model (default: 1)')
    prepare_parser.add_argument('--no-streme', action='store_true',
                                help='Stop after preparation (do not run STREME)')
    prepare_parser.add_argument('--nmotifs', type=int, default=200, help='STREME --nmotifs (default: 200)')
    prepare_parser.add_argument('--minw', type=int, default=6, help='STREME --minw (default: 6)')
    prepare_parser.add_argument('--maxw', type=int, default=20, help='STREME --maxw (default: 20)')
    prepare_parser.add_argument('--thresh', type=float, default=0.05, help='STREME --thresh (default: 0.05)')

    # Scan subcommand: FIMO every genome's STREME motifs against its promoters
    scan_parser = subparsers.add_parser(
        'scan',
        help='FIMO-scan every prepared genome\'s STREME motifs against its promoters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan all genomes in prepared/ at the default p-value threshold
  %(prog)s scan prepared/ --jobs 4

  # Use a q-value cutoff (recommended for typical analyses)
  %(prog)s scan prepared/ --thresh 0.05 --qv-thresh --jobs 4

  # Scan a different sequence set with each genome's STREME motifs
  %(prog)s scan prepared/ --sequence custom.fa --output custom_fimo/

For each streme_<genome>/ directory the scan picks: streme.txt as the motif
file, <genome>_prep/<genome>_promoters.masked as the default sequence file,
and <genome>_prep/background.txt as the background (else uses the motif
file's embedded background).
        """
    )
    scan_parser.add_argument('prepared_dir', help='Directory produced by `prepare` (contains streme_<genome>/ dirs)')
    scan_parser.add_argument('--output', '-o', help='Output root (default: same as prepared_dir)')
    scan_parser.add_argument('--jobs', '-j', type=int, default=1,
                             help='Genomes to scan in parallel (default: 1)')
    scan_parser.add_argument('--genomes', help='Comma-separated subset of genomes to scan')
    scan_parser.add_argument('--sequence', help='Override sequence FASTA for all genomes')
    scan_parser.add_argument('--bgfile', help='Override FIMO background for all genomes '
                                              '("motif-file" uses the motif file\'s embedded background)')
    scan_parser.add_argument('--thresh', type=float, default=1e-4,
                             help='FIMO match threshold (default: 1e-4)')
    scan_parser.add_argument('--qv-thresh', action='store_true',
                             help='Interpret --thresh as a q-value cutoff (e.g. 0.05)')
    scan_parser.add_argument('--no-qvalue', action='store_true', help='Skip FIMO q-value computation')
    scan_parser.add_argument('--max-strand', action='store_true',
                             help='Report only the higher-scoring strand of overlapping matches')
    scan_parser.add_argument('--parse-genomic-coord', action='store_true',
                             help='Parse genomic coords from FASTA headers')
    scan_parser.add_argument('--motif', help='Restrict to a specific motif ID')
    scan_parser.add_argument('--fimo-path', help='Path to the fimo executable')

    # Network subcommand: motif co-occurrence + gene clustering on a consolidated TSV
    network_parser = subparsers.add_parser(
        'network',
        help='Motif co-occurrence network + gene clustering on a consolidated TSV '
             '(exploratory analysis without expression data)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Defaults
  %(prog)s network outputs/consolidated_streme_sites.tsv --output network_results/

  # Tighter network + force 30 gene clusters
  %(prog)s network outputs/consolidated_streme_sites.tsv -o net/ \\
      --min-jaccard 0.2 --min-lift 3 --fdr 0.01 --n-clusters 30

Produces: motif_cooccurrence_edges.tsv, motif_cooccurrence_network.graphml
(import into Cytoscape/Gephi), motif_modules.tsv, gene_clusters.tsv,
gene_cluster_fingerprints.tsv, motif_positional_summary.tsv.
        """
    )
    network_parser.add_argument('consolidated_file', help='Consolidated STREME sites TSV')
    network_parser.add_argument('--output', '-o', default='network_results/',
                                help='Output directory (default: network_results/)')
    network_parser.add_argument('--min-motif-sites', type=int, default=10)
    network_parser.add_argument('--min-motifs-per-gene', type=int, default=2)
    network_parser.add_argument('--max-genes', type=int, default=5000)
    network_parser.add_argument('--min-jaccard', type=float, default=0.1)
    network_parser.add_argument('--min-lift', type=float, default=2.0)
    network_parser.add_argument('--fdr', type=float, default=0.05)
    network_parser.add_argument('--n-clusters', type=int)
    network_parser.add_argument('--cluster-distance', type=float, default=0.7)
    network_parser.add_argument('--skip-cooccurrence', action='store_true')
    network_parser.add_argument('--skip-clusters', action='store_true')
    network_parser.add_argument('--seed', type=int, default=42)

    # Network-viz subcommand: static figures from network outputs
    network_viz_parser = subparsers.add_parser(
        'network-viz',
        help='Create static PNG visualizations from motif-network outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default output: <network_dir>/figures/
    %(prog)s network_P_eatonii/

  # Custom figure output directory
    %(prog)s network_P_eatonii/ --output figures_P_eatonii/ --top-n 25
        """
    )
    network_viz_parser.add_argument('network_dir',
                                    help='Directory with motif_network.py outputs')
    network_viz_parser.add_argument('--output', '-o',
                                    help='Figure output directory '
                                         '(default: <network_dir>/figures/)')
    network_viz_parser.add_argument('--top-n', type=int, default=20,
                                    help='Top N categories in summary plots (default: 20)')
    network_viz_parser.add_argument('--network-nodes', type=int, default=120,
                                    help='Max nodes shown in network graphs (default: 120)')
    network_viz_parser.add_argument('--label-top', type=int, default=15,
                                    help='Label top N motifs in network graphs (default: 15)')
    network_viz_parser.add_argument('--dpi', type=int, default=200,
                                    help='Figure resolution in DPI (default: 200)')

    # Annotate subcommand: TOMTOM every prepared genome's STREME motifs vs a TF database
    annotate_parser = subparsers.add_parser(
        'annotate',
        help='TOMTOM-match every prepared genome\'s STREME motifs to a TF database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate vs PlantTFDB Arabidopsis (download separately, MEME format)
  %(prog)s annotate prepared/ --target-db ~/dbs/Ath_TF_binding_motifs.meme --jobs 4

  # E-value cutoff instead of the default q-value
  %(prog)s annotate prepared/ --target-db <db> --thresh 0.01 --evalue --jobs 4

Common databases (download in MEME format):
  - PlantTFDB 5.0 (Arabidopsis)
  - JASPAR plants (https://jaspar.genereg.net/)
  - CIS-BP (http://cisbp.ccbr.utoronto.ca/)
        """
    )
    annotate_parser.add_argument('prepared_dir', help='Directory produced by `prepare`')
    annotate_parser.add_argument('--target-db', required=True,
                                  help='Reference motif database (MEME format)')
    annotate_parser.add_argument('--output', '-o', help='Output root (default: same as prepared_dir)')
    annotate_parser.add_argument('--jobs', '-j', type=int, default=1,
                                  help='Genomes to annotate in parallel (default: 1)')
    annotate_parser.add_argument('--genomes', help='Comma-separated subset of genomes to annotate')
    annotate_parser.add_argument('--thresh', type=float, default=0.1,
                                  help='Significance threshold (q-value by default; '
                                       'E-value with --evalue) (default: 0.1)')
    annotate_parser.add_argument('--evalue', action='store_true',
                                  help='Interpret --thresh as an E-value cutoff')
    annotate_parser.add_argument('--no-ssc', action='store_true',
                                  help='Disable small-sample correction')
    annotate_parser.add_argument('--min-overlap', type=int, default=5,
                                  help='Minimum column overlap (default: 5)')
    annotate_parser.add_argument('--dist',
                                  choices=['pearson', 'ed', 'kullback', 'sandelin', 'allr'],
                                  default='pearson',
                                  help='Column-similarity metric (default: pearson)')
    annotate_parser.add_argument('--tomtom-path', help='Path to the tomtom executable')

    # Consolidate subcommand
    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate STREME motifs')
    consolidate_parser.add_argument('streme_dir', help='Directory with STREME results')
    consolidate_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    consolidate_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    consolidate_parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    # Validate subcommand
    validate_parser = subparsers.add_parser('validate', help='Validate motif consolidation quality')
    validate_parser.add_argument('consolidated_file', help='Consolidated motif TSV file (consolidated_streme_sites.tsv)')
    
    # Analyze subcommand
    analyze_parser = subparsers.add_parser('analyze', help='Run motif-expression analysis')
    analyze_parser.add_argument('consolidated_file', help='Consolidated motif TSV file')
    analyze_parser.add_argument('expression_file', help='Expression data file (long or wide format)')
    analyze_parser.add_argument('--type', choices=['absolute', 'relative'], default='relative',
                               help='Analysis type: absolute (per-genome) or relative (vs reference genome)')
    analyze_parser.add_argument('--reference-genome', '--reference-line', '-r',
                               dest='reference_line', default='IM767',
                               help='Reference genome for relative analysis (default: IM767). '
                                    "'--reference-line' is a deprecated alias.")
    analyze_parser.add_argument('--output', '-o', default='analysis_results/',
                               help='Output directory for analysis results')
    analyze_parser.add_argument('--detailed', action='store_true',
                               help='Use detailed motif features (not just presence/absence)')
    analyze_parser.add_argument('--top-motifs', type=int,
                               help='Only use top N most important motifs')
    analyze_parser.add_argument('--selection-method', choices=['frequency', 'variance', 'expression_corr'], 
                               default='variance',
                               help='Method for selecting top motifs: frequency (conservative), variance (differential), expression_corr (correlated)')
    
    # Gene-specific analysis subcommand
    gene_parser = subparsers.add_parser('gene-specific', help='Run gene-specific motif analysis')
    gene_parser.add_argument('consolidated_file', help='Consolidated motif TSV file')
    gene_parser.add_argument('expression_file', help='Expression data file')
    gene_parser.add_argument('--genes', help='Comma-separated list of target genes')
    gene_parser.add_argument('--gene-file', help='File containing one gene ID per line')
    gene_parser.add_argument('--type', choices=['absolute', 'relative'], default='relative',
                            help='Analysis type (default: relative)')
    gene_parser.add_argument('--reference-genome', '--reference-line', '-r',
                            dest='reference_line', default='IM767',
                            help='Reference genome for relative analysis (default: IM767). '
                                 "'--reference-line' is a deprecated alias.")
    gene_parser.add_argument('--output', '-o', default='gene_specific_results/',
                            help='Output directory')
    
    
    # Full pipeline subcommand
    full_parser = subparsers.add_parser('full', help='Run complete pipeline')
    full_parser.add_argument('streme_dir', help='Directory with STREME results '
                                                '(or where prepare writes them when --manifest is given)')
    full_parser.add_argument('expression_file', nargs='?', help='Expression data file for analysis step')
    full_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    full_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    full_parser.add_argument('--analysis-type', choices=['absolute', 'relative'], default='relative',
                            help='Analysis type for expression analysis step')
    full_parser.add_argument('--reference-genome', '--reference-line', '-r',
                            dest='reference_line', default='IM767',
                            help='Reference genome for relative analysis (default: IM767). '
                                 "'--reference-line' is a deprecated alias.")
    # Optional genome-preparation front stage (runs into streme_dir first).
    full_parser.add_argument('--manifest',
                            help='Run prepare first from this genome manifest, writing '
                                 'streme_<genome>/ dirs into streme_dir')
    full_parser.add_argument('--jobs', '-j', type=int, default=1,
                            help='Genomes to prepare in parallel when --manifest is used')
    full_parser.add_argument('--threads', type=int, default=1,
                            help='Threads per heavy step when --manifest is used')
    full_parser.add_argument('--upstream', '-u', type=int, default=1000,
                            help='Promoter upstream bp when --manifest is used (default: 1000)')
    full_parser.add_argument('--mask', choices=['none', 'repeatmasker', 'dust'], default='repeatmasker',
                            help='Masking step when --manifest is used (default: repeatmasker)')
    full_parser.add_argument('--species', help='Species for RepeatMasker when --manifest is used')
    full_parser.add_argument('--mask-lib',
                            help='Custom RepeatMasker library FASTA when --manifest is used '
                                 '(per-genome override via manifest `lib` column; takes precedence over --species)')
    full_parser.add_argument('--extra-mask', choices=['none', 'dust', 'trf', 'both'],
                            default='none',
                            help='Chain dust/TRF after the primary masker when --manifest is used')
    full_parser.add_argument('--trf-path', help='Path to the trf executable when --manifest is used')
    full_parser.add_argument('--chromosomes',
                            help='Restrict promoters to these sequences when --manifest is used '
                                 '(comma-separated names or a file; per-genome overrides via manifest)')
    full_parser.add_argument('--contig-pattern',
                            help='Restrict promoters to sequences matching this regex when --manifest is used')
    full_parser.add_argument('--masker-path',
                            help='Path to RepeatMasker/dust executable when --manifest is used')
    full_parser.add_argument('--fasta-get-markov-path',
                            help='Path to fasta-get-markov executable when --manifest is used')
    full_parser.add_argument('--streme-path',
                            help='Path to streme executable when --manifest is used')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Get script directory and project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    if args.command == 'prepare':
        success = run_prepare(args)
        sys.exit(0 if success else 1)

    elif args.command == 'scan':
        success = run_scan(args)
        sys.exit(0 if success else 1)

    elif args.command == 'annotate':
        success = run_annotate(args)
        sys.exit(0 if success else 1)

    elif args.command == 'network':
        cmd = [
            sys.executable, str(project_root / 'cli_tools' / 'motif_network.py'),
            args.consolidated_file,
            '--output', args.output,
            '--min-motif-sites', str(args.min_motif_sites),
            '--min-motifs-per-gene', str(args.min_motifs_per_gene),
            '--max-genes', str(args.max_genes),
            '--min-jaccard', str(args.min_jaccard),
            '--min-lift', str(args.min_lift),
            '--fdr', str(args.fdr),
            '--cluster-distance', str(args.cluster_distance),
            '--seed', str(args.seed),
        ]
        if args.n_clusters:
            cmd += ['--n-clusters', str(args.n_clusters)]
        if args.skip_cooccurrence:
            cmd.append('--skip-cooccurrence')
        if args.skip_clusters:
            cmd.append('--skip-clusters')
        success = run_command(cmd, "Running motif network + gene clustering")
        sys.exit(0 if success else 1)

    elif args.command == 'network-viz':
        cmd = [
            sys.executable, str(project_root / 'cli_tools' / 'motif_network_visualizer.py'),
            args.network_dir,
            '--top-n', str(args.top_n),
            '--network-nodes', str(args.network_nodes),
            '--label-top', str(args.label_top),
            '--dpi', str(args.dpi),
        ]
        if args.output:
            cmd += ['--output', args.output]
        success = run_command(cmd, "Creating network visualizations")
        sys.exit(0 if success else 1)

    elif args.command == 'consolidate':
        cmd = [
            'python', str(project_root / 'cli_tools' / 'streme_sites_consolidator.py'),
            'consolidate',
            args.streme_dir,
            '--output', args.output,
            '--threshold', str(args.threshold)
        ]
        if args.verbose:
            cmd.append('--verbose')
        
        success = run_command(cmd, "Consolidating STREME motifs")
        
    elif args.command == 'validate':
        cmd = [
            'python', str(project_root / 'cli_tools' / 'streme_sites_consolidator.py'),
            'validate',
            args.consolidated_file
        ]
        
        success = run_command(cmd, "Validating motif consolidation quality")
        
    elif args.command == 'analyze':
        # Choose the appropriate analyzer based on type
        if args.type == 'absolute':
            analyzer_script = 'motif_expression_analyzer.py'
        else:  # relative
            analyzer_script = 'relative_motif_analyzer.py'
        
        cmd = [
            'python', str(project_root / 'cli_tools' / analyzer_script),
            args.consolidated_file,
            args.expression_file,
            '--output', args.output
        ]
        
        if args.detailed:
            cmd.append('--detailed')
        if args.top_motifs:
            cmd.extend(['--top-motifs', str(args.top_motifs)])
        if args.type == 'relative':
            cmd.extend(['--reference-line', args.reference_line])
            cmd.extend(['--selection-method', args.selection_method])
        
        success = run_command(cmd, f"Running {args.type} motif-expression analysis")
        
    elif args.command == 'gene-specific':
        cmd = [
            'python', str(project_root / 'cli_tools' / 'gene_specific_analyzer.py'),
            args.consolidated_file,
            args.expression_file,
            '--type', args.type,
            '--output', args.output
        ]
        
        if args.genes:
            cmd.extend(['--genes', args.genes])
        if args.gene_file:
            cmd.extend(['--gene-file', args.gene_file])
        if args.type == 'relative':
            cmd.extend(['--reference-line', args.reference_line])
        
        success = run_command(cmd, f"Running gene-specific {args.type} analysis")
        
    elif args.command == 'full':
        # Optional Step 0: prepare genomes -> STREME results into streme_dir
        if args.manifest:
            prepare_args = argparse.Namespace(
                manifest=args.manifest, genome=None, fasta=None, annotation=None,
                expression=None, output=args.streme_dir, jobs=args.jobs, threads=args.threads,
                upstream=args.upstream, downstream=0, feature_type='gene',
                avoid_overlap=False, min_length=1,
                chromosomes=args.chromosomes, contig_pattern=args.contig_pattern,
                mask=args.mask, species=args.species, mask_lib=args.mask_lib,
                extra_mask=args.extra_mask, trf_path=args.trf_path,
                masker_path=args.masker_path,
                fasta_get_markov_path=args.fasta_get_markov_path,
                streme_path=args.streme_path,
                no_background=False, background_order=1, no_streme=False,
                nmotifs=200, minw=6, maxw=20, thresh=0.05,
                run_fimo=False, fimo_thresh=1e-4, fimo_qv_thresh=False,
                fimo_no_qvalue=False, fimo_max_strand=False, fimo_motif=None,
                fimo_path=None,
                run_tomtom=False, tomtom_db=None, tomtom_thresh=0.1,
                tomtom_evalue=False, tomtom_no_ssc=False,
                tomtom_min_overlap=5, tomtom_dist='pearson', tomtom_path=None,
            )
            if not run_prepare(prepare_args):
                print("⚠️  Some genomes failed to prepare; continuing with what succeeded.")

        # Step 1: Consolidate motifs
        cmd1 = [
            'python', str(project_root / 'cli_tools' / 'streme_sites_consolidator.py'),
            'consolidate',
            args.streme_dir,
            '--output', args.output,
            '--threshold', str(args.threshold)
        ]
        
        if not run_command(cmd1, "Step 1: Consolidating STREME motifs"):
            sys.exit(1)
        
        # Step 2: Validate consolidation
        consolidated_file = os.path.join(args.output, 'consolidated_streme_sites.tsv')
        if os.path.exists(consolidated_file):
            cmd2 = [
                'python', str(project_root / 'cli_tools' / 'streme_sites_consolidator.py'),
                'validate',
                consolidated_file
            ]
            if not run_command(cmd2, "Step 2: Validating consolidation quality"):
                print("⚠️  Validation failed, but continuing...")
        else:
            print("⚠️  Consolidated file not found, skipping validation")
        
        # Step 3: Expression analysis (if expression file provided)
        if args.expression_file and os.path.exists(consolidated_file):
            if args.analysis_type == 'absolute':
                analyzer_script = 'motif_expression_analyzer.py'
            else:  # relative
                analyzer_script = 'relative_motif_analyzer.py'
            
            analysis_output = os.path.join(args.output, f'{args.analysis_type}_analysis_results')
            cmd3 = [
                'python', str(project_root / 'cli_tools' / analyzer_script),
                consolidated_file,
                args.expression_file,
                '--output', analysis_output
            ]
            
            if args.analysis_type == 'relative':
                cmd3.extend(['--reference-line', args.reference_line])
            
            run_command(cmd3, f"Step 3: Running {args.analysis_type} motif-expression analysis")
            
            
            print(f"\n🎉 Full pipeline completed!")
            print(f"Consolidation results: {args.output}")
            print(f"Analysis results: {analysis_output}")
        else:
            if not args.expression_file:
                print(f"\n🎉 Consolidation and validation completed!")
                print(f"Results are in: {args.output}")
                print(f"\nTo run expression analysis:")
                print(f"python {sys.argv[0]} analyze {consolidated_file} expression.tsv --type relative")
            else:
                print("⚠️  Cannot run expression analysis - consolidated file missing")
        
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
