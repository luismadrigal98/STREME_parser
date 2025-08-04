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
  1. consolidate - Consolidate STREME motifs across lines
  2. map-genes   - Map motifs to individual genes
  3. full        - Run complete pipeline

Examples:
  # Step 1: Consolidate motifs
  %(prog)s consolidate /path/to/streme/results --output outputs/
  
  # Step 2: Map motifs to genes for specific line
  %(prog)s map-genes outputs/consolidated_motif_catalog.txt sequences.fasta IM502
  
  # Full pipeline
  %(prog)s full /path/to/streme/results /path/to/sequences/ --lines IM502,IM664,IM767
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline steps')
    
    # Consolidate subcommand
    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate STREME motifs')
    consolidate_parser.add_argument('streme_dir', help='Directory with STREME results')
    consolidate_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    consolidate_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    consolidate_parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    # Map genes subcommand
    map_parser = subparsers.add_parser('map-genes', help='Map motifs to genes')
    map_parser.add_argument('motif_catalog', help='Consolidated motif catalog')
    map_parser.add_argument('gene_sequences', help='Gene sequences FASTA file')
    map_parser.add_argument('line_name', help='Line name (e.g., IM502)')
    map_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    map_parser.add_argument('--detailed', action='store_true', help='Detailed report')
    
    # Full pipeline subcommand
    full_parser = subparsers.add_parser('full', help='Run complete pipeline')
    full_parser.add_argument('streme_dir', help='Directory with STREME results')
    full_parser.add_argument('sequences_dir', help='Directory with gene sequence files')
    full_parser.add_argument('--lines', required=True, help='Comma-separated list of line names')
    full_parser.add_argument('--output', '-o', default='outputs/', help='Output directory')
    full_parser.add_argument('--threshold', '-t', type=float, default=0.75, help='Similarity threshold')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Get script directory
    script_dir = Path(__file__).parent
    
    if args.command == 'consolidate':
        cmd = [
            'python', str(script_dir / 'cli_tools' / 'motif_consolidator.py'),
            args.streme_dir,
            '--output', args.output,
            '--threshold', str(args.threshold)
        ]
        if args.verbose:
            cmd.append('--verbose')
        
        success = run_command(cmd, "Consolidating STREME motifs")
        
    elif args.command == 'map-genes':
        cmd = [
            'python', str(script_dir / 'cli_tools' / 'gene_motif_mapper.py'),
            args.motif_catalog,
            args.gene_sequences,
            args.line_name,
            '--output', args.output
        ]
        if args.detailed:
            cmd.append('--detailed')
        
        success = run_command(cmd, f"Mapping motifs to genes for {args.line_name}")
        
    elif args.command == 'full':
        lines = args.lines.split(',')
        
        print(f"Running full pipeline for {len(lines)} lines: {lines}")
        
        # Step 1: Consolidate motifs
        cmd1 = [
            'python', str(script_dir / 'cli_tools' / 'motif_consolidator.py'),
            args.streme_dir,
            '--output', args.output,
            '--threshold', str(args.threshold)
        ]
        
        if not run_command(cmd1, "Step 1: Consolidating STREME motifs"):
            sys.exit(1)
        
        # Step 2: Map motifs to genes for each line
        catalog_file = os.path.join(args.output, 'consolidated_motif_catalog.txt')
        
        for line in lines:
            # Find sequence file for this line
            seq_file = None
            for ext in ['.fasta', '.fa', '.fna']:
                potential_file = os.path.join(args.sequences_dir, f"{line}{ext}")
                if os.path.exists(potential_file):
                    seq_file = potential_file
                    break
                potential_file = os.path.join(args.sequences_dir, f"Genes_{line}_DNA_lifted{ext}")
                if os.path.exists(potential_file):
                    seq_file = potential_file
                    break
            
            if not seq_file:
                print(f"⚠️  Could not find sequence file for line {line}")
                continue
            
            cmd2 = [
                'python', str(script_dir / 'cli_tools' / 'gene_motif_mapper.py'),
                catalog_file,
                seq_file,
                line,
                '--output', args.output,
                '--detailed'
            ]
            
            if not run_command(cmd2, f"Step 2: Mapping motifs to genes for {line}"):
                print(f"⚠️  Failed to process line {line}, continuing with others...")
                continue
        
        print(f"\n🎉 Full pipeline completed!")
        print(f"Results are in: {args.output}")
        
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
