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
    Optional column:  expression
    """
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
                "expression": (row[cols["expression"]].strip()
                               if "expression" in cols and row.get(cols["expression"]) else None),
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
        promoters = work_dir / f"{genome}_promoters.fasta"
        genome_prep.extract_promoters(
            spec["fasta"], spec["annotation"], str(promoters),
            upstream=opts["upstream"], downstream=opts["downstream"],
            feature_type=opts["feature_type"], avoid_overlap=opts["avoid_overlap"],
            min_length=opts["min_length"], genome_name=genome, verbose=True,
        )
        result["steps"]["promoters"] = str(promoters)
        streme_input = promoters

        if opts["mask"] != "none":
            masked = genome_prep.mask_sequences(
                str(streme_input),
                output=str(work_dir / f"{genome}_promoters.masked"),
                masker=opts["mask"], species=opts["species"], threads=opts["threads"],
            )
            result["steps"]["masked"] = masked
            streme_input = Path(masked)

        background = None
        if opts["background"]:
            background = genome_prep.build_background(
                str(streme_input), output=str(work_dir / "background.txt"),
                order=opts["background_order"],
            )
            result["steps"]["background"] = background

        if opts["run_streme"]:
            streme_dir = out_root / f"streme_{genome}"
            genome_prep.run_streme(
                str(streme_input), str(streme_dir),
                nmotifs=opts["nmotifs"], minw=opts["minw"], maxw=opts["maxw"],
                thresh=opts["thresh"], background=background, threads=opts["threads"],
            )
            result["steps"]["streme"] = str(streme_dir)
    except Exception as exc:  # noqa: BLE001 - surface per-genome failure to caller
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def run_prepare(args):
    """Drive genome preparation for one or many genomes, parallel across genomes."""
    if args.manifest:
        specs = read_genome_manifest(args.manifest)
    else:
        if not (args.genome and args.fasta and args.annotation):
            print("Error: provide --manifest OR all of --genome, --fasta, --annotation")
            return False
        specs = [{
            "genome": args.genome, "fasta": args.fasta,
            "annotation": args.annotation, "expression": args.expression,
        }]

    opts = {
        "output": args.output,
        "upstream": args.upstream, "downstream": args.downstream,
        "feature_type": args.feature_type, "avoid_overlap": args.avoid_overlap,
        "min_length": args.min_length,
        "mask": args.mask, "species": args.species,
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
    # Masking / background / STREME
    prepare_parser.add_argument('--mask', choices=['none', 'repeatmasker', 'dust'],
                                default='repeatmasker', help='Masking step (default: repeatmasker)')
    prepare_parser.add_argument('--species', help='Species for RepeatMasker')
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
                avoid_overlap=False, min_length=1, mask=args.mask, species=args.species,
                no_background=False, background_order=1, no_streme=False,
                nmotifs=200, minw=6, maxw=20, thresh=0.05,
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
