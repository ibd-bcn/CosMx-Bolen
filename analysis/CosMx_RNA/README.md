

# CosMx RNA Analysis – README

## General Information

Generated on **18 November 2025**  
Last modified on **18 November 2025**

This directory contains all scripts, intermediate data structures, and outputs associated with the **CosMx RNA** workflow used in the study *Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease*.  
CosMx RNA data were acquired in **2023** using NanoString CosMx™ Spatial Molecular Imaging.

The workflow was executed using:

- R version 4.3 or higher
- Seurat version 5
- InSituType
- Harmony
- ggplot2
- Python version 3.8 or higher (for cell–cell interaction)

Environment dependencies are provided at the repository root.

---

## File Overview

This directory contains the complete CosMx RNA analysis workflow.  
Scripts follow a numerical prefix that corresponds to the processing order.

### Files and Subdirectories

```
0_curation.R
1_annotation_supervised.R
2.normalization.R
3_cell_interaction.R
abundance_enrichment.R
scotia_cell_int.py
Files/
Markers/
Molecules/
Objects/
Polygons/
Results/
README.md
```


---

### File Naming Conventions

- `0_*` → Data loading, object construction, quality control  
- `1_*` → Supervised annotation  
- `2_*` → Normalization and downstream integration  
- `3_*` → Cell–cell interaction  
- Additional scripts → enrichment, abundance calculations, utilities

### File Formats Included

- CSV files (metadata, marker lists)
- RDS files (Seurat objects)
- TXT files (gene lists)
- PNG and PDF files (QC plots, marker visualizations)
- PY scripts (Python-based interaction analysis)

### Relationship Between Files

- `0_curation.R` creates the initial Seurat objects in `Objects/`
- `1_annotation_supervised.R` loads curated objects and applies supervised annotation
- `2.normalization.R` performs normalization and saves results in `Results/`
- `3_cell_interaction.R` and `scotia_cell_int.py` compute cell–cell interactions
- Intermediate QC files, polygons, molecules, and metadata are stored in their respective subdirectories

---

## Data-Specific Information

### 1. Raw data stored in /data/

Contains intermediate CosMx-derived tables used during preprocessing.

Typical files include:

- `exprMat.csv` – gene-by-cell expression matrix
- `metadata.csv` – cell-level metadata
- `tx_file.csv` – transcript-level molecule table

Common columns:

- `cell_ID` – unique cell identifier  
- `target_name` – gene target  
- `x_global_px`, `y_global_px` – global spatial coordinates  
- `Area` – segmentation polygon area  

Missing values are encoded as `NA`.

---

### 2. /analysis/CosMx_RNA/Molecules/

Contains molecule-level raw input exported from CosMx.

Typical columns:

- `TargetName` – gene name  
- `QScore` – detection confidence  
- `x_local`, `y_local` – cell-local coordinates  
- `x_global`, `y_global` – global coordinates  
- `CellID` – segmentation assignment  

---

### 3. /analysis/CosMx_RNA/Polygons/

Contains segmentation polygons for each cell.  
Used for visualization, spatial QC, and interaction analysis.

---

### 4. analysis/CosMx_RNA/Molecules/

Includes all Seurat objects generated during the workflow:

- `qc_seurat.RDS` – unfiltered raw object  
- `qc_seurat_pass.RDS` – object after QC filtering  
- `annotated_seurat.RDS` – object with supervised cell type labels  
- `normalized_seurat.RDS` – final normalized object for analysis

Assays include:

- `RNA` – gene expression  
- `NegProbe` – negative probe counts  
- `SystemControl` – system control probes  

---

### 5. analysis/CosMx_RNA/Markers/

Contains marker gene lists and marker-based QC plots.  
Files include `markers_global.csv` and `*_markers.png`.

---

### 6. analysis/CosMx_RNA/Results/

Contains all exported results used for figure generation:

- Abundance summaries  
- Enrichment matrices  
- Normalized counts  
- Cell-level annotation tables  

### 7. analysis/CosMx_RNA/Files/

Contains intermediate files. Irrelevant. 

---

## Methodology

### Instruments and Data Collection

Data were produced using **NanoString CosMx™ Spatial RNA** (2023).  
Raw spatial tables were exported from CosMx nCounter software.

### Software and Tools

- Seurat (preprocessing, QC, annotation)
- InSituType (supervised annotation)
- Harmony (integration)
- ggplot2 (visualization)
- Python scotia (cell–cell interaction)

---

## Workflow Summary

### Step 0 – Curation and QC (`0_curation.R`)
- Load raw tables  
- Create Seurat object  
- Create RNA, Negative Probe, and System Control assays  
- Apply QC filters:

  - low counts  
  - high negative probe proportion  
  - low gene complexity  
  - polygon area filtering (Grubbs test)  
  - target-level QC  
- Save QC-filtered objects to `Objects/`

---

### Step 1 – Supervised Annotation (`1_annotation_supervised.R`)
- Load scRNA-seq reference  
- Select markers  
- Run InSituType supervised classification  
- Export annotated objects and probability tables  

---

### Step 2 – Normalization (`2.normalization.R`)
- Apply SCTransform or log-normalization  
- Perform scaling and PCA  
- Apply Harmony integration (if used)  
- Export normalized objects to `Results/`

---

### Step 3 – Cell–Cell Interaction  
Files:  
- `3_cell_interaction.R`  
- `scotia_cell_int.py`

Steps:

1. Export coordinates and cell types  
2. Run Python-based scotia interaction analysis  
3. Import results into R for summarization and visualization  

---

### Step 4 – Abundance and Enrichment  
File: `abundance_enrichment.R`

Computes:

- differential abundance  
- enrichment scores  
- summary matrices  

---

## Data Access and Sharing

### Publications Using This Dataset

This dataset supports the analysis described in:  
*Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease.*

### Recommended Citation

Bolen et al., 2025. *Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease.*

### License

MIT License.

### Funding Acknowledgment

“This research was funded in whole or in part by Aligning Science Across Parkinson’s (ASAP Grant #) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).”

