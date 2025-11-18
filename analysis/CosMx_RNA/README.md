# CosMx RNA Analysis – README

## General Information

Generated on **18 November 2025**  
Last modified on **18 November 2025**

This directory contains all scripts, intermediate data structures, and outputs associated with the **CosMx RNA** workflow used in the study *Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease*.

CosMx RNA data were acquired in **2023** using NanoString CosMx™ Spatial Molecular Imaging.

---

## File Overview

This directory contains the full CosMx RNA analysis workflow.  
Scripts follow a numerical prefix corresponding to the order of processing.

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

## File Naming Conventions

- `0_*` → data loading, QC, object creation  
- `1_*` → supervised annotation  
- `2_*` → normalization and integration  
- `3_*` → cell–cell interaction  
- `abundance_enrichment.R` → abundance and enrichment summaries  

### File Formats

- `.csv` → metadata, molecule tables, marker lists  
- `.RDS` → Seurat objects  
- `.png`, `.pdf` → visualizations  
- `.py` → Python scripts (scotia interaction analysis)  

---

## Data-Specific Information

### Raw Data (`/data/`)

CosMx-exported RNA tables per FOV:

- `exprMat.csv` – gene-by-cell expression matrix  
- `metadata.csv` – cell-level metadata  
- `tx_file.csv` – transcript-level molecule table  
- `polygons.csv` – cell boundaries  

Important metadata columns:

- `cell`  
- `fov`  
- `Area`  
- `x_global_px`, `y_global_px`  
- `tissue`, `patient`

Negative probes:

- `Negative*`  
- `SystemControl*`

Missing values are encoded as `NA`.

---

## Directory: `Objects/`

Contains all Seurat objects:

- `qc_seurat.RDS` – unfiltered object  
- `qc_seurat_pass.RDS` – QC-filtered object  
- `annotated_seurat.RDS` – object with supervised annotation  
- `normalized_seurat.RDS` – final normalized object  

Assays included:

- **RNA** – gene expression  
- **NegProbe** – negative probe counts  
- **SystemControl** – system control probes  

---

## Directory: `Molecules/`

Contains molecule-level raw input exported from CosMx.

Columns:

- `TargetName` – gene name  
- `QScore` – detection confidence  
- `x_local`, `y_local` – cell-local coordinates  
- `x_global`, `y_global` – global coordinates  
- `CellID` – segmentation assignment  

---

## Directory: `Polygons/`

One CSV per FOV containing:

- Cell IDs  
- Polygon vertex coordinates  
- Cell boundaries for visualization and spatial analysis


---

## Directory: `Markers/`

Contains marker gene lists and marker-based QC plots:

- `markers_global.csv`  
- `*_markers.png`  

---

## Directory: `Results/`

Contains final outputs:

- UMAP and PCA embeddings  
- Cell type annotation tables  
- Abundance and enrichment matrices  
- Cell–cell interaction results  
- Visualizations for manuscript figures

---

## Directory: `Files/`

Intermediate files not used for downstream analysis.

---

## Workflow

### Step 0 – Curation & QC (`0_curation.R`)

- Load raw RNA tables  
- Build Seurat object (`RNA`, `NegProbe`, `SystemControl` assays)  
- Append metadata (cell ID, FOV)  
- Apply QC filters:

  **QC Flag 1 – Low counts**  
  **QC Flag 2 – High negative probe proportion**  
  **QC Flag 3 – Low gene complexity**  
  **QC Flag 4 – Polygon area outliers (Grubbs test)**  
  **QC Flag 5 – Target-level QC**  

Outputs saved to `/Objects/`.

---

### Step 1 – Supervised Annotation (`1_annotation_supervised.R`)

- Load scRNA-seq reference  
- Select marker genes  
- Run InSituType supervised classification  
- Export annotated objects and probability tables

---

### Step 2 – Normalization (`2.normalization.R`)

- Apply SCTransform or log-normalization  
- Perform scaling and PCA  
- Optional Harmony integration  
- Export normalized objects to `/Results/`

---

### Step 3 – Cell–Cell Interaction

Files:

- `3_cell_interaction.R`  
- `scotia_cell_int.py`

Steps:

1. Export coordinates and cell types  
2. Run Python-based scotia interaction analysis  
3. Import results into R for summarization and visualization  
4. Save outputs to `/Results/`

---

### Step 4 – Abundance & Enrichment (`abundance_enrichment.R`)

Computes:

- Cell type abundance per tissue and patient  
- Enrichment matrices  
- Condition-specific differences  

Outputs stored in `/Results/`.

---

## Citation

If using this dataset or code, please cite:
Bolen et al., 2025.  
*Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease.*

---

## Funding

“This research was funded by Aligning Science Across Parkinson’s (ASAP Grant #) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).”

---

## License

MIT License.
