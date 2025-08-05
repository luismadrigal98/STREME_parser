# STREME Sites Consolidator - Output Columns Guide

This document explains each column in the consolidated TSV file produced by the STREME Sites Consolidator tool.

## Output File Columns

The consolidated TSV file contains the following 16 columns, providing comprehensive information about motif occurrences across genetic lines:

### 1. `consolidated_motif_id`
**Format:** `MOTIF_001`, `MOTIF_002`, etc.
**Description:** A unique identifier assigned to each consolidated motif cluster. Motifs with similar consensus patterns (ignoring numeric prefixes from STREME) are grouped together and assigned the same consolidated ID. This allows tracking of the same regulatory element across different genetic lines.

### 2. `original_motif_id`
**Format:** `1-AAAAWMWTTTT`, `2-CGCGCGCG`, etc.
**Description:** The original motif identifier from the STREME `sites.tsv` file. Includes the numeric prefix (e.g., "1-", "2-") that STREME assigns to distinguish motifs within a single analysis. This preserves traceability to the original STREME output.

### 3. `motif_consensus`
**Format:** `AAAAWMWTTTT`, `CGCGCGCG`, etc.
**Description:** The consensus sequence pattern for this motif, with the numeric prefix removed. This is the actual DNA sequence pattern that defines the motif, using IUPAC nucleotide codes (W=A/T, M=A/C, etc.). This is what's used for clustering similar motifs.

### 4. `line`
**Format:** `DPR_14`, `SWB_2`, `YOD_7`, etc.
**Description:** The genetic line identifier, extracted from the directory name containing the STREME analysis. Represents different Mimulus guttatus lines or populations being analyzed.

### 5. `gene_id`
**Format:** `Migut.A00001.v2.0`, `Migut.B01234.v2.0`, etc.
**Description:** The unique identifier for the gene containing this motif occurrence. Corresponds to gene annotations in the Mimulus guttatus genome assembly.

### 6. `start_pos`
**Format:** Integer (e.g., `156`, `2341`)
**Description:** The starting position of the motif occurrence within the gene sequence (1-based indexing). Indicates where the regulatory element begins relative to the gene start.

### 7. `end_pos`
**Format:** Integer (e.g., `167`, `2352`)
**Description:** The ending position of the motif occurrence within the gene sequence (1-based indexing). Together with `start_pos`, defines the exact location of the regulatory element.

### 8. `strand`
**Format:** `+` or `-`
**Description:** The DNA strand on which the motif was found. `+` indicates the forward (sense) strand, `-` indicates the reverse (antisense) strand. Important for understanding regulatory directionality.

### 9. `score`
**Format:** Decimal number (e.g., `8.45`, `12.3`)
**Description:** The STREME score for this motif occurrence, indicating the strength or confidence of the match. Higher scores generally indicate better matches to the motif consensus pattern.

### 10. `sequence`
**Format:** `AAAATCATTTT`, `CGCGCGCG`, etc.
**Description:** The actual DNA sequence found at this position that matches the motif pattern. This is the real sequence from the genome, which may contain slight variations from the consensus pattern.

### 11. `representative_sequence`
**Format:** `AAAAWMWTTTT`, `CGCGCGCG`, etc.
**Description:** The consensus pattern of the most common motif variant within this consolidated cluster. Used as the "canonical" representation of this motif type across all lines.

### 12. `relative_position`
**Format:** Integer (e.g., `1`, `2`, `3`)
**Description:** The sequential position of this motif within the gene (per line), ordered by genomic position. The first motif in a gene has `relative_position = 1`, the second has `relative_position = 2`, etc. Useful for analyzing regulatory element organization.

### 13. `total_motifs_in_gene`
**Format:** Integer (e.g., `3`, `7`, `12`)
**Description:** The total number of motif occurrences found within this gene in this specific line. Indicates the regulatory complexity or density for this gene.

### 14. `relative_position_fraction`
**Format:** Decimal (e.g., `0.33`, `0.67`, `1.0`)
**Description:** The relative position as a fraction of total motifs in the gene (`relative_position / total_motifs_in_gene`). Values range from 0 to 1, where 0.33 means this is the first motif out of 3 total motifs in the gene. Useful for comparing motif positions across genes with different motif counts.

### 15. `cluster_size`
**Format:** Integer (e.g., `2`, `5`, `8`)
**Description:** The number of original STREME motif variants that were consolidated into this `consolidated_motif_id`. A higher number indicates this motif pattern was found with more sequence variations across different lines or within the same analysis.

### 16. `length`
**Format:** Integer (e.g., `8`, `11`, `15`)
**Description:** The length of the motif sequence in base pairs. Calculated as `end_pos - start_pos + 1`. Indicates the size of the regulatory element.

## Usage Examples

### Finding genes with multiple regulatory elements:
```bash
# Genes with more than 5 motifs in any line
awk -F'\t' '$13 > 5' consolidated_motifs.tsv | cut -f4,5 | sort | uniq
```

### Analyzing motif distribution patterns:
```bash
# Most common motif types
cut -f1 consolidated_motifs.tsv | sort | uniq -c | sort -nr
```

### Comparing regulatory complexity across lines:
```bash
# Average motifs per gene by line
awk -F'\t' 'NR>1 {sum[$4] += $13; count[$4]++} END {for(line in sum) print line, sum[line]/count[line]}' consolidated_motifs.tsv
```

## Quality Control Notes

- **Missing Data:** Empty cells indicate data was not available in the original STREME output
- **Consistency:** All positions and coordinates are 1-based to match biological conventions
- **Sorting:** Output is sorted by consolidated_motif_id, then line, then gene_id, then start_pos for easy analysis
- **Validation:** Use `cluster_size` and `total_motifs_in_gene` to identify potentially interesting regulatory regions

## Related Files

- `consolidated_motifs_summary.txt` - Statistical summary of the consolidation results
- Original STREME `sites.tsv` files in respective line directories
- `requirements.txt` - Python dependencies (pandas)
