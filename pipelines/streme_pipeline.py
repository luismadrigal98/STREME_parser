#!/usr/bin/env python3
"""
STREME Analysis Pipeline - Main Pipeline Steps:
  1. consolidate      - Consolidate STREME motifs across lines
  2. validate         - Validate motif consolidation quality
  3. analyze          - Run motif-expression analysis (absolute or relative)
  4. full             - Run complete pipeline (steps 1-3)

This orchestrates the complete pipeline from STREME output consolidation
through validation to expression analysis.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

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
  4. full            - Run complete pipeline (steps 1-3)

Examples:
  # Step 1: Consolidate motifs
  %(prog)s consolidate /path/to/streme/results --output outputs/
  
  # Step 2: Validate consolidation quality
  %(prog)s validate outputs/consolidated_streme_sites.tsv
  
  # Step 3: Expression analysis
  %(prog)s analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --output results/
  
  # Full pipeline
  %(prog)s full /path/to/streme/results expression.tsv --output outputs/ --analysis-type relative

Expression Analysis Types:
  - absolute: Per-line analysis (each line analyzed independently)  
  - relative: Comparative analysis (each line compared to IM767 baseline)
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline steps')
    
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
                               help='Analysis type: absolute (per-line) or relative (vs IM767)')
    analyze_parser.add_argument('--output', '-o', default='analysis_results/',
                               help='Output directory for analysis results')
    analyze_parser.add_argument('--detailed', action='store_true',
                               help='Use detailed motif features (not just presence/absence)')
    analyze_parser.add_argument('--top-motifs', type=int,
                               help='Only use top N most important motifs')
    
    # Full pipeline subcommand
    full_parser = subparsers.add_parser('full', help='Run complete pipeline')
    full_parser.add_argument('streme_dir', help='Directory with STREME results')
    full_parser.add_argument('expression_file', nargs='?', help='Expression data file for analysis step')
    full_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    full_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    full_parser.add_argument('--analysis-type', choices=['absolute', 'relative'], default='relative',
                            help='Analysis type for expression analysis step')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Get script directory and project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    if args.command == 'consolidate':
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
        
        success = run_command(cmd, f"Running {args.type} motif-expression analysis")
        
    elif args.command == 'full':
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
