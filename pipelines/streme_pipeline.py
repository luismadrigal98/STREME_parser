#!/usr/bin/env python3
"""
STREME Analysis Pipeline

Comprehensive toolkit for analyzing STREME motif discovery results and their
effects on gene expression across multiple genetic lines.

Commands:
- consolidate: Consolidate STREME motifs across lines
- validate: Validate motif consolidation quality
- extract-features: Extract ML features from motifs
- analyze-expression: Basic motif-expression analysis
- analyze-mixed: Mixed-effects analysis (RECOMMENDED)
- analyze-comprehensive: Advanced regulatory analysis
- full: Run complete pipeline

Statistical Approach:
Uses mixed-effects models (the gold standard in regulatory genomics) to:
- Handle gene-specific baseline differences (random intercepts)
- Allow motif effects to vary by gene (random slopes)
- Borrow statistical strength across genes
- Provide population-level conclusions about motif effects

For detailed help: streme-parser <command> --help
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n{'='*60}")
    print(f"🚀 {description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"❌ {description} failed: {e}")
        return False

def get_script_paths():
    """Get paths to CLI tools"""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    cli_tools = project_root / 'cli_tools'
    return project_root, cli_tools

def consolidate_command(args):
    """Run motif consolidation"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'motif_consolidator.py'),
        args.streme_dir,
        '--output', args.output,
        '--threshold', str(args.threshold)
    ]
    if args.verbose:
        cmd.append('--verbose')
    
    return run_command(cmd, "Consolidating STREME motifs across lines")

def validate_command(args):
    """Run motif validation"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'validate_consolidation.py'),
        args.consolidated_file
    ]
    if args.output:
        cmd.extend(['--output', args.output])
    
    return run_command(cmd, "Validating motif consolidation quality")

def extract_features_command(args):
    """Run feature extraction"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'motif_to_regression_features.py'),
        args.motif_file,
        args.fasta_file,
        '--output', args.output
    ]
    
    if args.simple:
        cmd.append('--simple')
    if args.top_motifs:
        cmd.extend(['--top-motifs', str(args.top_motifs)])
    if args.min_sites:
        cmd.extend(['--min-sites', str(args.min_sites)])
    
    return run_command(cmd, "Extracting motif features for analysis")

def analyze_expression_command(args):
    """Run basic expression analysis"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'motif_expression_analyzer.py'),
        args.motif_features,
        args.expression_data,
        '--output', args.output
    ]
    
    if args.method:
        cmd.extend(['--method', args.method])
    if args.top_motifs:
        cmd.extend(['--top-motifs', str(args.top_motifs)])
    
    return run_command(cmd, "Analyzing motif effects on gene expression")

def analyze_mixed_command(args):
    """Run mixed-effects analysis (RECOMMENDED approach)"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'mixed_effects_analyzer.py'),
        args.motif_features,
        args.expression_data,
        '--approach', args.approach,
        '--output', args.output
    ]
    
    if args.min_genes:
        cmd.extend(['--min-genes', str(args.min_genes)])
    
    return run_command(cmd, "Mixed-Effects Analysis (Gold Standard)")

def analyze_comprehensive_command(args):
    """Run comprehensive regulatory analysis"""
    project_root, cli_tools = get_script_paths()
    cmd = [
        sys.executable,
        str(cli_tools / 'comprehensive_motif_analyzer.py'),
        args.motif_features,
        args.expression_data,
        '--output', args.output
    ]
    
    if args.effects:
        cmd.extend(['--effects'] + args.effects)
    if args.position_window:
        cmd.extend(['--position-window', str(args.position_window)])
    
    return run_command(cmd, "Comprehensive Regulatory Analysis")

def full_pipeline_command(args):
    """Run the complete pipeline"""
    print("\n" + "="*80)
    print("🔬 STREME COMPLETE ANALYSIS PIPELINE")
    print("="*80)
    print("Steps: Consolidate → Validate → Extract Features → Mixed-Effects Analysis")
    print("="*80)
    
    success_count = 0
    total_steps = 0
    
    # Step 1: Consolidate motifs (if input directory provided)
    if args.input_dir:
        total_steps += 1
        consolidate_args = argparse.Namespace(
            streme_dir=args.input_dir,
            output=args.output_dir,
            threshold=args.similarity,
            verbose=False
        )
        if consolidate_command(consolidate_args):
            success_count += 1
            # Set motif file path for next steps
            args.motif_file = os.path.join(args.output_dir, 'consolidated_streme_sites.tsv')
        else:
            print("❌ Pipeline failed at consolidation step")
            return False
    
    # Step 2: Validate consolidation
    if args.motif_file and os.path.exists(args.motif_file):
        total_steps += 1
        validate_args = argparse.Namespace(
            consolidated_file=args.motif_file,
            output=args.output_dir
        )
        if validate_command(validate_args):
            success_count += 1
    
    # Step 3: Extract features
    if args.motif_file and args.fasta_file:
        total_steps += 1
        features_output = os.path.join(args.output_dir, 'motif_features.tsv')
        extract_args = argparse.Namespace(
            motif_file=args.motif_file,
            fasta_file=args.fasta_file,
            output=features_output,
            simple=args.simple,
            top_motifs=args.top_motifs,
            min_sites=None
        )
        if extract_features_command(extract_args):
            success_count += 1
            args.features_file = features_output
        else:
            print("❌ Pipeline failed at feature extraction step")
            return False
    
    # Step 4: Mixed-effects analysis (RECOMMENDED)
    if hasattr(args, 'features_file') and args.expression_data:
        total_steps += 1
        analysis_output = os.path.join(args.output_dir, 'mixed_effects_analysis')
        mixed_args = argparse.Namespace(
            motif_features=args.features_file,
            expression_data=args.expression_data,
            approach='all',  # Compare all approaches
            output=analysis_output,
            min_genes=50
        )
        if analyze_mixed_command(mixed_args):
            success_count += 1
    
    # Pipeline summary
    print(f"\n{'='*80}")
    print(f"📊 PIPELINE SUMMARY")
    print(f"{'='*80}")
    print(f"Steps completed successfully: {success_count}/{total_steps}")
    
    if success_count == total_steps:
        print("🎉 Pipeline completed successfully!")
        print(f"\n📁 Results location: {args.output_dir}")
        print("📈 Key outputs:")
        if args.motif_file:
            print(f"  • Consolidated motifs: {args.motif_file}")
        if hasattr(args, 'features_file'):
            print(f"  • Feature matrix: {args.features_file}")
        if total_steps >= 4:
            print(f"  • Statistical analysis: {analysis_output}/")
        print("\n💡 Next steps:")
        print("  • Review the mixed_effects_motif_coefficients.tsv for significant motifs")
        print("  • Check analysis_comparison.png for model performance")
        print("  • See docs/STATISTICAL_MODELING_GUIDE.md for interpretation help")
    else:
        print("⚠️  Pipeline completed with some errors.")
        print("Check the output above for specific error messages.")
    
    print(f"{'='*80}")
    return success_count == total_steps

def main():
    parser = argparse.ArgumentParser(
        description='STREME Analysis Pipeline - Comprehensive motif-expression analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available Commands:
  consolidate         Consolidate STREME motifs across lines
  validate            Validate motif consolidation quality
  extract-features    Extract features for analysis
  analyze-expression  Basic motif-expression analysis
  analyze-mixed       Mixed-effects analysis (RECOMMENDED)
  analyze-comprehensive  Advanced regulatory analysis
  full                Complete pipeline

Statistical Approaches:
  • Mixed-effects models (analyze-mixed): Gold standard in regulatory genomics
  • Gene-by-gene analysis: Individual gene modeling
  • Hierarchical analysis: Groups similar genes

Examples:
  # Complete pipeline (recommended)
  %(prog)s full --input-dir streme_results/ --fasta sequences.fa --expression expr.tsv --output results/
  
  # Mixed-effects analysis only (if you have features already)
  %(prog)s analyze-mixed motif_features.tsv expression.tsv --output analysis/
  
  # Individual steps
  %(prog)s consolidate streme_results/ --output consolidated/
  %(prog)s extract-features consolidated_motifs.tsv sequences.fa --simple
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline commands')
    
    # Consolidate command
    consolidate_parser = subparsers.add_parser(
        'consolidate', 
        help='Consolidate STREME motifs across genetic lines'
    )
    consolidate_parser.add_argument('streme_dir', help='Directory containing STREME results')
    consolidate_parser.add_argument('--output', '-o', default='outputs/', 
                                   help='Output directory')
    consolidate_parser.add_argument('--threshold', '-t', type=float, default=0.8,
                                   help='Motif similarity threshold (0-1)')
    consolidate_parser.add_argument('--verbose', '-v', action='store_true',
                                   help='Verbose output')
    
    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate motif consolidation quality'
    )
    validate_parser.add_argument('consolidated_file', 
                                help='Consolidated motif file (consolidated_streme_sites.tsv)')
    validate_parser.add_argument('--output', '-o', help='Output directory')
    
    # Extract features command
    features_parser = subparsers.add_parser(
        'extract-features',
        help='Extract features from consolidated motifs'
    )
    features_parser.add_argument('motif_file', help='Consolidated motif file')
    features_parser.add_argument('fasta_file', help='FASTA sequence file')
    features_parser.add_argument('--output', '-o', default='motif_features.tsv',
                                help='Output feature file')
    features_parser.add_argument('--simple', action='store_true',
                                help='Extract only presence/absence features')
    features_parser.add_argument('--top-motifs', type=int,
                                help='Use only top N most frequent motifs')
    features_parser.add_argument('--min-sites', type=int, default=10,
                                help='Minimum sites required for motif inclusion')
    
    # Basic expression analysis command
    expression_parser = subparsers.add_parser(
        'analyze-expression',
        help='Basic motif-expression analysis'
    )
    expression_parser.add_argument('motif_features', help='Motif feature file')
    expression_parser.add_argument('expression_data', help='Expression data file')
    expression_parser.add_argument('--output', '-o', default='expression_analysis',
                                  help='Output directory')
    expression_parser.add_argument('--method', choices=['linear', 'random_forest'],
                                  default='linear', help='Analysis method')
    expression_parser.add_argument('--top-motifs', type=int,
                                  help='Use only top N motifs')
    
    # Mixed-effects analysis command (RECOMMENDED)
    mixed_parser = subparsers.add_parser(
        'analyze-mixed',
        help='Mixed-effects analysis (RECOMMENDED approach)'
    )
    mixed_parser.add_argument('motif_features', help='Motif feature file')
    mixed_parser.add_argument('expression_data', help='Expression data file')
    mixed_parser.add_argument('--approach', 
                             choices=['mixed', 'gene-by-gene', 'hierarchical', 'all'],
                             default='all', help='Statistical approach')
    mixed_parser.add_argument('--output', '-o', default='mixed_effects_results',
                             help='Output directory')
    mixed_parser.add_argument('--min-genes', type=int, default=50,
                             help='Minimum genes required for mixed effects')
    
    # Comprehensive analysis command
    comprehensive_parser = subparsers.add_parser(
        'analyze-comprehensive',
        help='Comprehensive regulatory analysis'
    )
    comprehensive_parser.add_argument('motif_features', help='Motif feature file')
    comprehensive_parser.add_argument('expression_data', help='Expression data file')
    comprehensive_parser.add_argument('--output', '-o', default='comprehensive_analysis',
                                     help='Output directory')
    comprehensive_parser.add_argument('--effects', nargs='+',
                                     choices=['presence', 'position', 'variation'],
                                     default=['presence', 'position', 'variation'],
                                     help='Effects to analyze')
    comprehensive_parser.add_argument('--position-window', type=int, default=100,
                                     help='Window size for position analysis')
    
    # Full pipeline command
    full_parser = subparsers.add_parser(
        'full',
        help='Run complete analysis pipeline'
    )
    full_parser.add_argument('--input-dir', help='STREME results directory')
    full_parser.add_argument('--motif-file', help='Pre-consolidated motif file')
    full_parser.add_argument('--fasta-file', required=True, help='FASTA sequence file')
    full_parser.add_argument('--expression-data', required=True, help='Expression data file')
    full_parser.add_argument('--output-dir', '-o', default='pipeline_results',
                            help='Output directory for all results')
    full_parser.add_argument('--similarity', type=float, default=0.8,
                            help='Motif similarity threshold')
    full_parser.add_argument('--simple', action='store_true',
                            help='Use simple features (presence/absence only)')
    full_parser.add_argument('--top-motifs', type=int,
                            help='Use only top N motifs')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Route to appropriate command function
    success = False
    
    if args.command == 'consolidate':
        success = consolidate_command(args)
    elif args.command == 'validate':
        success = validate_command(args)
    elif args.command == 'extract-features':
        success = extract_features_command(args)
    elif args.command == 'analyze-expression':
        success = analyze_expression_command(args)
    elif args.command == 'analyze-mixed':
        success = analyze_mixed_command(args)
    elif args.command == 'analyze-comprehensive':
        success = analyze_comprehensive_command(args)
    elif args.command == 'full':
        success = full_pipeline_command(args)
    else:
        parser.print_help()
        return 1
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())