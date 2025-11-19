# CosMx Protein Analysis – README

## General Information

Generated on **18 November 2025**  
Last modified on **19 November 2025**

This directory contains all scripts, intermediate data structures, and outputs associated with the **CosMx Protein** workflow used in the study *Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease*.

CosMx Protein data were acquired in **2023–2024** using NanoString CosMx™ Spatial Molecular Imaging with the **68-plex protein panel**.

---

## File Overview

This directory contains the full CosMx Protein analysis workflow.  
Scripts follow a numerical prefix corresponding to the order of processing.

```
0_curation.R
1_normalization.R
2.0_celltyping.R
2.1_celltyping.py
2.2_celltyping.R
3_cell_interaction.R
abundances_enrichment.R
scotia_cell_int.py
Objects/
Polygons/
Files/
Results/
README.md
```

---

## File Naming Conventions

- `0_*` → data loading, QC, object creation  
- `1_*` → normalization  
- `2_*` → cross-modal cell type annotation (RNA → Protein)  
- `3_*` → spatial cell–cell interaction  
- `abundances_enrichment.R` → abundance and enrichment summaries  

### File Formats

- `.csv` → metadata, correspondence lists, counts  
- `.RDS` → Seurat objects  
- `.png`, `.pdf` → visualizations  
- `.py` → Python scripts (scotia interaction analysis)

---

## Data-Specific Information

### Raw Data (`/data/`)

CosMx-exported protein tables per FOV:

- `exprMat_file.csv` – protein intensity matrix  
- `metadata_file.csv` – segmentation metadata  
- `polygons.csv` – cell boundaries  

Important metadata columns:

- `cell`  
- `fov`  
- `Area`  
- `x_global_px`, `y_global_px`  
- `tissue`, `patient`

Negative probes:

- `Rb IgG`  
- `Ms IgG1`

Missing values are encoded as `NA`.

---

## Directory: `Objects/`

Contains all Seurat objects:

- `seurats.RDS`  
- `qc_seurats.RDS`  
- `qc_seurats_filt.RDS`  
- `norm.RDS`  
- `sc_ref_cut.RDS`  
- `counts.RDS`  
- `meta.csv`  
- `counts_SC_cut.csv`, `meta_SC_cut.csv` (MaxFuse inputs)

Assays included:

- **Prot** – protein counts  
- **Negprob** – negative control probe intensities

---

## Directory: `Polygons/`

One CSV per FOV containing:

- Cell IDs  
- Polygon vertex coordinates  
- Cell boundaries for visualization and spatial network analysis

---

## Directory: `Results/`

Contains final outputs:

- UMAP and PCA embeddings  
- Cell type annotation tables  
- Abundance and enrichment tables  
- Cell–cell interaction results  
- Visualizations for manuscript figures

---

## Directory: `Files/`

Intermediate files not used for downstream analysis.

---

## Workflow

### Step 0 – Curation & QC (`0_curation.R`)

- Load raw protein tables  
- Build Seurat object (`Prot` + `Negprob` assays)  
- Append metadata (patient, tissue, FOV)  
- Export polygons  
- Generate QC flags:

  **QC Flag 1 – Protein count distribution**  
  **QC Flag 2 – Negative probes** 
  **QC Flag 3 – Polygon area**  

Outputs saved to `Objects/`.

---

### Step 1 – Normalization (`1_normalization.R`)

Procedure:

1. Total intensity normalization  
2. Arcsinh transform (cofactor = 50)  
3. Variable feature selection  
4. Scaling  
5. PCA (40 PCs)  
6. UMAP (PCs 1–25)

Output: `norm.RDS`

---

### Step 2 – Cell Typing

#### Step 2.0 – Reference Setup (`2.0_celltyping.R`)

- Load Trigos scRNA-seq reference  
- Filter relevant compartments  
- Downsample to 2,000 cells per label  
- Build protein-to-gene correspondence table  
- Export matrices for MaxFuse

#### Step 2.1 – MaxFuse (`2.1_celltyping.py`)

- Perform RNA ↔ Protein cross-modal mapping  
- Generate label transfer file

#### Step 2.2 – Annotation Cleanup (`2.2_celltyping.R`)

- Import MaxFuse predictions  
- Apply final annotation  
- Save updated object

---

### Step 3 – Cell–Cell Interactions

Files:

- `3_cell_interaction.R`  
- `scotia_cell_int.py`

Workflow:

1. Export cell coordinates + cell types  
2. Run **scotia** (Python)  
3. Import edge list and interaction scores  
4. Summaries by tissue, patient, and condition  
5. Save outputs to `/Results/`

---

### Step 4 – Abundance & Enrichment (`abundances_enrichment.R`)

Computes:

- Cell type abundance per tissue and patient  
- Enrichment matrices  
- Condition-specific differences  

All saved under `/Results/`.

---

## Citation

If using this dataset or code, please cite:

Bolen et al., 2025.  
*Spatial single-cell multiomics reveals peripheral immune dysfunction in Parkinson’s and inflammatory bowel disease.*

---

## Funding

“This research was funded by Aligning Science Across Parkinson’s (ASAP-020527) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).”

---

## License

MIT License.
