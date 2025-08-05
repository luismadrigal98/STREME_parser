#!/usr/bin/env python3
"""
MEME Analysis Pipeline - Global CLI Tool

This tool provides a unified interface for all MEME/STREME analysis tools.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

def get_script_dir():
    """Get the directory containing this script"""
    return Path(__file__).parent.absolute()

def run_consolidator(args):
    """Run the STREME sites consolidator"""
    script_path = get_script_dir() / "cli_tools" / "streme_sites_consolidator.py"
    
    if not script_path.exists():
        print(f"Error: Consolidator script not found at {script_path}")
        return 1
    
    # Build command
    cmd = [sys.executable, str(script_path)]
    
    # Add positional argument
    cmd.append(args.streme_results_dir)
    
    # Add optional arguments
    if args.output:
        cmd.extend(['--output', args.output])
    if args.threshold:
        cmd.extend(['--threshold', str(args.threshold)])
    if args.lines:
        cmd.extend(['--lines', args.lines])
    if args.overlap_threshold:
        cmd.extend(['--overlap-threshold', str(args.overlap_threshold)])
    if args.no_merge_overlaps:
        cmd.append('--no-merge-overlaps')
    if args.verbose:
        cmd.append('--verbose')
    
    # Run the command
    print(f"Running: {' '.join(cmd)}")
    return subprocess.call(cmd)

def run_motif_mapper(args):
    """Run the gene motif mapper (placeholder for future implementation)"""
    print("Gene motif mapper - Coming soon!")
    print("This will map consolidated motifs to genes and create regulatory maps")
    return 0

def run_sequence_grouper(args):
    """Run sequence grouping by motifs (placeholder for future implementation)"""
    print("Sequence grouper - Coming soon!")
    print("This will group sequences by their motif patterns")
    return 0

def main():
    """Main CLI interface"""
    parser = argparse.ArgumentParser(
        description='MEME Analysis Pipeline - Unified CLI for motif analysis tools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available Commands:
  consolidate    Consolidate STREME sites.tsv files across multiple lines
  map-genes      Map consolidated motifs to genes (coming soon)
  group-seqs     Group sequences by motif patterns (coming soon)

Examples:
  # Consolidate STREME results with overlap merging
  %(prog)s consolidate /path/to/streme/results/ --output comprehensive_analysis
  
  # Consolidate specific lines with custom thresholds
  %(prog)s consolidate /path/to/streme/results/ \\
    --lines IM502,IM664 \\
    --threshold 0.8 \\
    --overlap-threshold 0.7 \\
    --output targeted_analysis
  
  # Disable overlap merging for debugging
  %(prog)s consolidate /path/to/streme/results/ \\
    --no-merge-overlaps \\
    --output debug_analysis

Directory Structure Expected:
  streme_results/
  ├── streme_IM502/sites.tsv
  ├── streme_IM664/sites.tsv
  └── streme_IM767/sites.tsv
        """
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='COMMAND'
    )
    
    # Consolidator subcommand
    consolidate_parser = subparsers.add_parser(
        'consolidate',
        help='Consolidate STREME sites.tsv files across multiple lines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Consolidate STREME sites.tsv files and create comprehensive motif maps',
        epilog="""
This tool processes STREME sites.tsv files from multiple genetic lines,
consolidates similar motifs, merges overlapping hits, and creates
comprehensive gene-level regulatory maps.

Key Features:
- Consolidates similar motifs using IUPAC-aware comparison
- Merges overlapping motif hits (removes STREME sliding window artifacts)
- Handles repetitive sequences (poly-A, poly-T) appropriately
- Provides exact genomic coordinates and sequences
- Generates comprehensive statistics and summaries
        """
    )
    
    consolidate_parser.add_argument(
        'streme_results_dir',
        help='Directory containing STREME output folders with sites.tsv files'
    )
    
    consolidate_parser.add_argument(
        '--output', '-o',
        help='Output file prefix (default: consolidated_streme_sites)'
    )
    
    consolidate_parser.add_argument(
        '--threshold', '-t',
        type=float,
        help='Similarity threshold for motif consolidation (default: 0.75)'
    )
    
    consolidate_parser.add_argument(
        '--lines', '-l',
        help='Comma-separated list of specific lines to process (e.g., IM502,IM664)'
    )
    
    consolidate_parser.add_argument(
        '--overlap-threshold',
        type=float,
        help='Minimum overlap fraction to merge motifs (default: 0.5)'
    )
    
    consolidate_parser.add_argument(
        '--no-merge-overlaps',
        action='store_true',
        help='Disable merging of overlapping motif hits'
    )
    
    consolidate_parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    # Gene mapper subcommand (placeholder)
    map_parser = subparsers.add_parser(
        'map-genes',
        help='Map consolidated motifs to genes (coming soon)',
        description='Create gene-level regulatory maps from consolidated motifs'
    )
    
    map_parser.add_argument(
        'consolidated_file',
        help='Consolidated motifs TSV file'
    )
    
    map_parser.add_argument(
        '--output', '-o',
        help='Output file prefix'
    )
    
    # Sequence grouper subcommand (placeholder)
    group_parser = subparsers.add_parser(
        'group-seqs',
        help='Group sequences by motif patterns (coming soon)',
        description='Group sequences based on their regulatory motif patterns'
    )
    
    group_parser.add_argument(
        'consolidated_file',
        help='Consolidated motifs TSV file'
    )
    
    group_parser.add_argument(
        '--output', '-o',
        help='Output directory'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if command was provided
    if not args.command:
        parser.print_help()
        return 1
    
    # Route to appropriate function
    if args.command == 'consolidate':
        return run_consolidator(args)
    elif args.command == 'map-genes':
        return run_motif_mapper(args)
    elif args.command == 'group-seqs':
        return run_sequence_grouper(args)
    else:
        print(f"Unknown command: {args.command}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
