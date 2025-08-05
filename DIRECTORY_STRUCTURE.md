# Directory Structure

This directory has been organized into the following structure:

## Root Level

- `README.md` - Main project documentation
- `DIRECTORY_STRUCTURE.md` - This file
- `main.py` - Main entry point for the entire pipeline
- `requirements.txt` - Python package dependencies

## Organized Directories

### `/scripts/`

Contains executable scripts and utilities:

- `MEME_val.py` - MEME validation script
- `consolidate_streme_results.sh` - Shell script for consolidating STREME results
- `remote_streme_all_lines.sh` - Remote STREME execution script
- `run_analysis.sh` - Main analysis runner script

### `/docs/`

Contains documentation:

- `MOTIF_ANALYSIS_GUIDE.md` - Comprehensive guide for motif analysis

### `/pipelines/`

Contains pipeline definitions and main pipeline scripts:

- `meme_pipeline.py` - Main MEME pipeline implementation
- `Pipeline_promotor_discovery.yaml` - Pipeline configuration file

### `/cli_tools/`

Contains command-line interface tools:

- `gene_motif_mapper.py` - Maps motifs to individual genes with position scoring
- `motif_consolidator.py` - Consolidates motif data

### `/archive/`

Contains older versions of files and test utilities:

- Historical analysis scripts kept for reference
- Test files: `test_motif_parsing.py`, `test_similarity.py`, `test_proper_similarity.py`
- Superseded Python scripts from development iterations

### `/outputs/`

Directory for storing analysis outputs:

- `.gitkeep` - Ensures directory is tracked in git

## Migration Notes

The directory was cleaned by:
1. Removing duplicate files that existed in both root and archive directories
2. Organizing files by function into logical subdirectories
3. Keeping the most recent versions of tools in their appropriate locations
4. Preserving all historical files in the archive directory

This structure makes it easier to:
- Find specific tools and scripts
- Understand the project organization
- Maintain and develop new features
- Run tests and documentation builds
