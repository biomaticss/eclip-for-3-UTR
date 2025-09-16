# eCLIP Data Analysis Workflow

> **Project Goal**: Analyze the binding sites of RNA-binding proteins (RBPs) `PTBP1`, `DDX3X`, `ELAVL1`, and `PABPN1` on genes `STAT5B`, `GSPT1`, `CCND2`, and `CTNNBIP1`, focusing on signal intensity near proximal polyadenylation sites (PAS) (±200 bp), and generate clear visualization plots.

---

## 1. Background and Objectives

This project utilizes eCLIP data from ENCODE to investigate the binding patterns of four RBPs on target genes, emphasizing regions near proximal PAS. The main objectives are:
- **Generate Intersection Files**: Intersect RBP BED files with PAS regions (`pas_nearby.bed`) to identify binding sites.
- **Extract Signal**: Extract signal intensity around PAS (±500 bp) from bigWig files.
- **Visualize**: Generate individual plots for each gene-RBP pair (e.g., `CCND2-PTBP1_binding.png`), displaying peak positions (scatter plot) and signal intensity (line plot).

The output consists of 16 plots (4 genes × 4 RBPs), saved in the `plots/` directory, designed for research handover and display on a homepage.

---

## 2. Data Preparation

### 2.1 Input Files
| File Name | Description | Source |
|-----------|-------------|--------|
| `pas_nearby.bed` | Proximal PAS coordinates (7 columns: chrom, start, end, name, score, strand, gene) | Generated from `human.PAS.txt` |
| `gene_coords.bed` | Gene coordinates (6 columns: chrom, start, end, name, score, strand) | Generated from `position.tsv` canonical transcripts |
| `PTBP1_K562.bed`, etc. | RBP peak data (BED format) | ENCODE |
| `PTBP1_K562_plus.bigWig`, etc. | RBP signal data (bigWig format, plus/minus strands) | ENCODE |
| `hg38.chrom.sizes` | Chromosome sizes | UCSC Genome Browser (hg38) |

#### 2.1.1 Generating `pas_nearby.bed`
- **Source**: The file is derived from `human.PAS.txt`, a dataset containing polyadenylation site annotations (likely from PolyASite or similar databases).
- **Process**:
  1. **Filter Proximal PAS**: A script (e.g., `select_proximal_pas.py`) was used to select the most proximal PAS for each gene in `human.PAS.txt`. For each gene (`STAT5B`, `GSPT1`, `CCND2`, `CTNNBIP1`), the PAS closest to the 3' end of the canonical transcript was chosen.
  2. **Format Conversion**: The selected PAS coordinates were formatted into BED format with 7 columns:
     - `chrom`: Chromosome (e.g., `chr17`).
     - `start`, `end`: PAS region coordinates (typically ±200 bp around the PAS).
     - `name`: PAS identifier (e.g., `chr17:40353638:-`).
     - `score`: Signal strength or confidence score (e.g., `2.59892720629`).
     - `strand`: Strand orientation (`+` or `-`).
     - `gene`: Associated gene name (e.g., `STAT5B`).
  3. **Output**: Saved as `pas_nearby.bed`.

**Example**:
```
chr17	42201419	42201820	chr17:40353638:-	2.59892720629	-	STAT5B
chr16	11872527	11872928	chr16:11966585:-	1.81640147122	-	GSPT1
chr12	4304552	4304953	chr12:4413919:+	2.09426857209	+	CCND2
chr1	9850282	9850683	chr1:9910541:-	2.12246050094	-	CTNNBIP1
```

#### 2.1.2 Generating `gene_coords.bed`
- **Source**: The file is derived from `position.tsv`, a table containing genomic coordinates of canonical transcripts for the target genes.
- **Process**:
  1. **Select Canonical Transcripts**: For each gene (`STAT5B`, `GSPT1`, `CCND2`, `CTNNBIP1`), the canonical transcript was identified from `position.tsv` (likely sourced from Ensembl or GENCODE).
  2. **Extract Coordinates**: Extracted chromosome, start, end, gene name, score (set to 0), and strand for each gene.
  3. **Format Conversion**: Converted to BED format with 6 columns:
     - `chrom`: Chromosome (e.g., `chr12`).
     - `start`, `end`: Gene boundaries (transcription start to end).
     - `name`: Gene name (e.g., `CCND2`).
     - `score`: Placeholder (set to 0).
     - `strand`: Strand orientation (`+` or `-`).
  4. **Output**: Saved as `gene_coords.bed`.

**Example**:
```
chr12	4273761	4305353	CCND2	0	+
chr1	9850449	9893237	CTNNBIP1	0	-
chr16	11873118	11916082	GSPT1	0	-
chr17	42199176	42276391	STAT5B	0	-
```

### 2.2 Environment Setup
Install dependencies:
```bash
pip install pandas matplotlib seaborn pybedtools pyBigWig
```

Ensure the `bigWigAverageOverBed` executable (from UCSC tools) is available.

---

## 3. Issues and Solutions

### 3.1 `pas_nearby.bed` Column Misalignment
**Issue**:
- `cat -A pas_nearby.bed` revealed extra spaces between `strand` and `gene` (e.g., `I-   STAT5B$`), causing pandas to parse it as 8 columns, with `gene` as `NaN`.
- The code's `pas_df['gene'].str.contains(gene, case=False, na=False)` failed to match gene names.

**Solution**:
Cleaned `pas_nearby.bed` to remove extra spaces and convert gene names to uppercase:
```bash
awk '{gsub(/[ \t]+/, "\t", $0); $7 = toupper($7); print}' OFS="\t" pas_nearby.bed > pas_nearby_clean.bed
mv pas_nearby_clean.bed pas_nearby.bed
```

Updated the code to use regex separator `sep=r'\t+'` and enforce `gene` as string type:
```python
pas_df = pd.read_csv(PAS_NEARBY_BED, sep=r'\t+', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'gene'],
                     dtype={'gene': str}, engine='python')
```

**Verification**:
```bash
cat -A pas_nearby.bed | head
cut -f7 pas_nearby.bed | sort | uniq
```
Output confirmed a 7-column structure and correct gene names (`CCND2`, `CTNNBIP1`, `GSPT1`, `STAT5B`).

### 3.2 Visualization Requirement Change
**Issue**:
- The original code generated summary plots (heatmap, bar chart, binding profile) combining all RBPs for each gene.
- New requirement: Generate individual plots for each gene-RBP pair (e.g., `CCND2-PTBP1_binding.png`), showing peak positions (scatter) and signal intensity (line).

**Solution**:
- Removed `plot_heatmap` and `plot_bar_chart` functions.
- Added `plot_single_gene_rbp` function to create individual plots, each focusing on PAS ±500 bp, including:
  - Scatter plot: Blue dots for peak positions.
  - Line plot: Green solid line (primary strand signal), orange dashed line (secondary strand signal).
  - Red dashed line: Proximal PAS position.

---
## 4. Final Code

The final code, modified from the original, generates 16 individual plots (4 genes × 4 RBPs), saved in the `plots/` directory.

eclip.py

---

## 5. Execution Steps

1. **Clean `pas_nearby.bed`**:
   ```bash
   awk '{gsub(/[ \t]+/, "\t", $0); $7 = toupper($7); print}' OFS="\t" pas_nearby.bed > pas_nearby_clean.bed
   mv pas_nearby_clean.bed pas_nearby.bed
   ```

2. **Verify Input Files**:
   ```bash
   ls PTBP1_K562.bed DDX3X_HepG2.bed ELAVL1_K562.bed PABPN1_HepG2.bed
   ls PTBP1_K562_plus.bigWig PTBP1_K562_minus.bigWig DDX3X_HepG2_plus.bigWig DDX3X_HepG2_minus.bigWig ELAVL1_K562_plus.bigWig ELAVL1_K562_minus.bigWig PABPN1_HepG2_plus.bigWig PABPN1_HepG2_minus.bigWig
   cat gene_coords.bed
   ```

3. **Run the Script**:
   ```bash
   python3 eclip.py
   ```

4. **Check Outputs**:
   - **Plots**: 16 images in `plots/` (e.g., `CCND2-PTBP1_binding.png`).
   - **Debug Info**: Confirm `Matches for {gene}` and `Generating plot for {gene}-{rbp}` outputs.
   - **Intersection Files**: e.g., `PTBP1_intersections.bed`.

---

## 6. Output Results

### 6.1 Plot Description
- **File Name**: `{gene}-{rbp}_binding.png` (e.g., `CCND2-PTBP1_binding.png`).
- **Content**:
  - **Scatter Plot**: Blue dots for peak positions (y=1).
  - **Line Plot**:
    - Green solid line: Primary strand signal (plus for + strand, minus for - strand, y=1.5 offset).
    - Orange dashed line: Secondary strand signal (minus for + strand, plus for - strand, y=2.0 offset).
  - **Red Dashed Line**: Proximal PAS position.
- **X-Axis**: PAS ±200 bp.
- **Y-Axis**: Peak positions and signal intensity.

### 6.2 Example Plot
(Embed example plot in your homepage)
```markdown
![CCND2-PTBP1 Binding](plots/CCND2-PTBP1_binding.png)
```

---

## 7. Notes

- **PAS Coordinate Validation**:
  - The PAS in `pas_nearby.bed` (e.g., `chr17:40353638:-`) may not lie within `STAT5B` (`chr17:42199176-42276391`) 3' UTR. Verify:
    ```bash
    grep "STAT5B" pas_nearby.bed
    ```
    Check in UCSC Genome Browser (hg38). If incorrect, revisit `human.PAS.txt` or the PAS selection script.

- **bigWig Files**:
  - Ensure all bigWig files are valid and use `chr*` naming (e.g., `chr1`, not `1`).
  - Verify:
    ```bash
    file ELAVL1_K562_plus.bigWig
    ```

- **Dependencies**:
  - Ensure `pyBigWig` is installed, or signal extraction will fail.
  - Confirm `bigWigAverageOverBed` is executable.

- **Debugging**:
  - If no plots are generated, check for `No PAS found` or `Warning: bigWig file not found` in output.
  - Verify `pas_df['gene']`:
    ```bash
    cut -f7 pas_nearby.bed | sort | uniq
    ```

---

## 8. Summary

This workflow successfully visualizes RBP binding sites and signal intensity, producing clear individual plots for research handover and homepage display. Future extensions could include:
- Adding bar charts to compare average signal intensity near PAS.
- Refining PAS selection to ensure 3' UTR localization.
- Incorporating additional RBPs or genes.

**Contact**: For further discussion, reach out via [hyut_2000@163.com].

---

*Generated on September 16, 2025*
