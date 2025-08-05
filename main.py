#!/usr/bin/env python3
"""
MEME Analysis Pipeline - Main Entry Point

This is the primary script for running the complete MEME motif analysis pipeline.
It provides a unified interface to all pipeline components.

Usage:
    python main.py --help
    python main.py consolidate /path/to/streme/results
    python main.py map-genes /path/to/catalog.txt /path/to/sequences/
"""

import os
import sys
import subprocess
from pathlib import Path

# Add the project directories to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root / "cli_tools"))
sys.path.insert(0, str(project_root / "pipelines"))

def main():
    """Main entry point - delegate to meme_pipeline.py"""
    pipeline_script = project_root / "pipelines" / "meme_pipeline.py"
    
    if not pipeline_script.exists():
        print(f"Error: Pipeline script not found at {pipeline_script}")
        sys.exit(1)
    
    # Pass all arguments to the pipeline script
    cmd = [sys.executable, str(pipeline_script)] + sys.argv[1:]
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user")
        sys.exit(1)

if __name__ == "__main__":
    main()
