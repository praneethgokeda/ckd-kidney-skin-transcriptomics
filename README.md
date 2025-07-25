# Shared Inflammatory Signatures Between Kidney and Skin in Chronic Kidney Disease

**Author**: Praneeth Gokeda  
**Type**: Systems Biology / RNA-Seq Project  
**Tools**: FastQC, Trim Galore, HISAT2, HTSeq, edgeR, STRING, Enrichr, Cytoscape

---

## Abstract

This project explores the shared inflammatory mechanisms in kidney and skin tissues in the context of Chronic Kidney Disease (CKD). RNA-Seq data analysis revealed differential gene expression and pathway enrichment, identifying **arachidonic acid metabolism** as a central pathway linking both tissues through inflammation-related genes such as **CYP2C**, **ALOX12**, and **PTGS**.

---

## 📁 Folder Structure

systems_project/
├── data/ → Gene count files (.count.txt)
├── results/ → Plots, PPI networks, KEGG images, enrichment CSVs
├── report/ → Final PDF report with figures and references


---

## Data Source

- 3 kidney + 3 skin RNA-Seq samples  
- Source: Human Protein Atlas, Uppsala Biobank  
- Original study: [PMID: 25613900](https://pubmed.ncbi.nlm.nih.gov/25613900/)

Raw FASTQ/BAM files are not included due to size. Preprocessing and mapping were performed as described in the report.

---

## Pipeline Overview

1. **Preprocessing**: FastQC + Trim Galore  
2. **Mapping**: HISAT2 (GRCh38 genome)  
3. **Counting**: HTSeq-count  
4. **DEG Analysis**: edgeR in R  
5. **Enrichment**: Enrichr (GO, KEGG), STRING + Cytoscape  
6. **Visualization**: PCA, MDS, BCV, heatmaps, PPI, pathway maps

---

## Key Results

- 10,043 DEGs identified
  - 5,373 upregulated in kidney
  - 4,670 upregulated in skin
- Common pathway: **Arachidonic acid metabolism**
- Shared inflammatory genes: *ALOX12, PTGS, CYP2C19 (skin)* and *CYP2C, CYP4A11 (kidney)*
- Psoriasis gene signature overlaps with CKD (adj. p = 4.047e-7)

---

## Report

View full findings and references in the PDF:  
[📘 Final Report](./report/Systems_project_final_report.pdf)

---

## Notes

- Raw data not uploaded. See data source above.
- This project is focused on downstream interpretation and visualization.
- R scripts not included.

---

## License

This work is shared for academic and educational purposes under the MIT License.