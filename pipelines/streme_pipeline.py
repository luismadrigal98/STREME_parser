#!/usr/bin/env python3
"""
STREME Analysis Pipeline - Main Pipeline Steps:
  1. consolidate      - Consolidate STREME motifs across lines
  2. validate         - Validate motif consolidation quality
  3. full             - Run complete pipeline (steps 1-2)

This orchestrates the complete pipeline from STREME output consolidation
to validation. For feature extraction, use the CLI tools directly.
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
        description='STREME Analysis Pipeline - Complete motif discovery and validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline Steps:
  1. consolidate     - Consolidate STREME motifs across lines
  2. validate        - Validate motif consolidation quality
  3. full            - Run complete pipeline (steps 1-2)

Examples:
  # Step 1: Consolidate motifs
  %(prog)s consolidate /path/to/streme/results --output outputs/
  
  # Step 2: Validate consolidation quality
  %(prog)s validate outputs/consolidated_streme_sites.tsv
  
  # Full pipeline
  %(prog)s full /path/to/streme/results --output outputs/ --threshold 0.75

For feature extraction, use the CLI tools directly:
  python cli_tools/streme_sites_consolidator.py features consolidated_streme_sites.tsv
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
            run_command(cmd2, "Step 2: Validating consolidation quality")
        else:
            print("⚠️  Consolidated file not found, skipping validation")
        
        print(f"\n🎉 Full pipeline completed!")
        print(f"Results are in: {args.output}")
        print(f"\nNext step: Extract features using:")
        print(f"python cli_tools/streme_sites_consolidator.py features {consolidated_file}")
        
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
