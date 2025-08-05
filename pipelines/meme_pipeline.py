#!/usr/bin/env python3
"""
MEME Analysis Pipeline - Master CLI tool for motif discovery and analysis.

This tool orchestrates the complete pipeline from STREME output consolidation
to gene-level motif mapping for expression analysis.
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
        description='MEME Analysis Pipeline - Complete motif discovery and mapping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline Steps:
  1. consolidate     - Consolidate STREME motifs across lines
  2. validate        - Validate motif consolidation quality
  3. extract-features - Extract regression features from consolidated motifs
  4. full            - Run complete pipeline (steps 1-3)

Examples:
  # Step 1: Consolidate motifs
  %(prog)s consolidate /path/to/streme/results --output outputs/
  
  # Step 2: Validate consolidation quality
  %(prog)s validate outputs/consolidated_streme_sites.tsv
  
  # Step 3: Extract features (detailed)
  %(prog)s extract-features outputs/consolidated_streme_sites.tsv --expression expr_data.tsv
  
  # Step 3: Extract features (simple presence/absence only)
  %(prog)s extract-features outputs/consolidated_streme_sites.tsv --simple
  
  # Full pipeline
  %(prog)s full /path/to/streme/results --output outputs/
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
    
    # Full pipeline subcommand
    full_parser = subparsers.add_parser('full', help='Run complete pipeline')
    full_parser.add_argument('streme_dir', help='Directory with STREME results')
    full_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    full_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    
    # Feature extraction subcommand
    features_parser = subparsers.add_parser('extract-features', help='Extract features from consolidated motif data')
    features_parser.add_argument('motif_file', help='Consolidated motif TSV file')
    features_parser.add_argument('--expression', '-e', help='Expression data file (TSV with line, gene, expression columns)')
    features_parser.add_argument('--top-motifs', '-t', type=int, help='Include only top N most frequent motifs')
    features_parser.add_argument('--min-sites', '-m', type=int, default=10, help='Minimum sites required for motif inclusion')
    features_parser.add_argument('--output-prefix', '-o', default='motif_regression', help='Output file prefix')
    features_parser.add_argument('--simple', '-s', action='store_true', help='Generate only presence/absence features (binary)')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Get script directory and project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    if args.command == 'consolidate':
        cmd = [
            'python', str(project_root / 'cli_tools' / 'motif_consolidator.py'),
            args.streme_dir,
            '--output', args.output,
            '--threshold', str(args.threshold)
        ]
        if args.verbose:
            cmd.append('--verbose')
        
        success = run_command(cmd, "Consolidating STREME motifs")
        
    elif args.command == 'validate':
        cmd = [
            'python3', str(project_root / 'validate_consolidation.py'),
            args.consolidated_file
        ]
        
        success = run_command(cmd, "Validating motif consolidation quality")
        
    elif args.command == 'full':
        # Step 1: Consolidate motifs
        cmd1 = [
            'python', str(project_root / 'cli_tools' / 'motif_consolidator.py'),
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
                'python3', str(project_root / 'validate_consolidation.py'),
                consolidated_file
            ]
            run_command(cmd2, "Step 2: Validating consolidation quality")
        else:
            print("⚠️  Consolidated file not found, skipping validation")
        
        print(f"\n🎉 Full pipeline completed!")
        print(f"Results are in: {args.output}")
        
    elif args.command == 'extract-features':
        cmd = [
            'python', str(project_root / 'cli_tools' / 'motif_to_regression_features.py'),
            args.motif_file,
            '--output-prefix', args.output_prefix,
            '--min-sites', str(args.min_sites)
        ]
        
        if args.expression:
            cmd.extend(['--expression', args.expression])
        if args.top_motifs:
            cmd.extend(['--top-motifs', str(args.top_motifs)])
        if args.simple:
            cmd.append('--simple')
        
        success = run_command(cmd, "Extracting regression features from motif data")
        
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
