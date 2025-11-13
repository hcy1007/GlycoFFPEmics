# Glyco FFPEmics Glycopeptide Analysis System

## System Overview

Glyco FFPEmics is a comprehensive platform designed for glycopeptide omics analysis of FFPE (formalin‑fixed paraffin‑embedded) samples. The system integrates a complete workflow from data upload to results visualization, supporting three core modules: glycan analysis, glycopeptide qualitative analysis, and glycopeptide quantitative analysis, providing powerful tools for glycopeptide research.

## System Features

- **Multi‑dimensional analysis**: Supports glycan analysis, glycopeptide qualitative analysis, and glycopeptide quantitative analysis.
- **Rich visualization**: Provides heatmaps, volcano plots, ROC curves, Venn diagrams, and more.
- **Functional enrichment**: Integrates GO and KEGG pathway enrichment analysis.
- **Machine learning**: Supports PCA, LDA, SVM, and other ML methods.
- **User‑friendly**: Intuitive web interface with support for downloading and batch exporting results.

## Application Scenarios

- Glycopeptide omics research on FFPE samples  
- Discovery of cancer glycosylation biomarkers  
- Differential analysis of glycopeptide quantification  
- Functional study of glycosylation pathways  

## System Architecture

### Core Modules

1. **Data Upload & Quality Control**
   - Supports uploading multiple data file formats  
   - Data quality check and sample matching validation  
   - Missing value analysis and data distribution visualization  

2. **Glycan Analysis**
   - Molecular weight distribution  
   - Identification of significantly different glycans  
   - Volcano plot visualization  

3. **Glycopeptide Qualitative Analysis**
   - Sample distribution analysis  
   - Correlation analysis  
   - PCA dimensionality reduction  
   - Glycosylation site distribution  
   - Sankey diagram analysis  
   - Unique statistics  

4. **Glycopeptide Quantitative Analysis (Core Module)**
   - Heatmap analysis (six display modes)  
   - Regional heatmap analysis  
   - Venn diagram analysis  
   - LDA plots and GO analysis  
   - KEGG analysis  
   - ROC analysis  
   - Volcano plot analysis  

## Installation Guide

### System Requirements

- **Operating System**: Windows 10/11, macOS 10.14+, Linux (Ubuntu 18.04+)  
- **Memory**: Recommended 8GB+  
- **Storage**: At least 2GB free space  
- **Network**: Stable internet connection (for enrichment analysis)  
- **Browser**: Chrome, Firefox, Safari, Edge (HTML5 supported)  

### Installation Steps

1. **Install R environment**
   ```bash
   # Download and install R (https://cran.r-project.org/)
   # Recommended version: R 4.2.0 or higher
   ```

2. **Install RStudio** (optional but recommended)
   ```bash
   # Download and install RStudio (https://posit.co/download/rstudio-desktop/)
   ```

3. **Clone or download the project**
   ```bash
   git clone [project URL]
   # Or download and unzip the ZIP package
   ```

4. **Install dependencies**
   ```r
   # Run in R
   source("install_dependencies.R")
   ```

5. **Launch the application**
   ```r
   # Run in R
   source("appchange.R")
   ```

## User Guide

### Quick Start

1. **Prepare data**
   - Group‑normalized data file (.csv)  
   - Patient information file (.csv)  
   - Glycan type file (.csv)  

2. **Start the system**
   - Run `appchange.R`  
   - Open the displayed URL in your browser  

3. **Upload data**
   - Upload the three required files on the “Data Upload” page  
   - Verify sample matching  

4. **Quality control**
   - Review data quality on the QC page  

5. **Select analysis module**
   - Choose the desired analysis module  
   - Run and view results  

### Detailed Analysis Workflow

#### Glycan Analysis
1. Upload data file  
2. Set valid value proportion threshold  
3. Run glycan analysis  
4. View molecular weight distribution and differential results  

#### Glycopeptide Qualitative Analysis
1. Upload data file  
2. Run qualitative analysis  
3. Review sample distribution, correlation, PCA, etc.  

#### Glycopeptide Quantitative Analysis
1. Upload specialized quantitative files  
2. Select analysis submodule  
3. Run and view results  

## Data Format Requirements

### Required Files

1. **Group‑normalized data file (.csv)**
   - Contains glycopeptide identification and quantification data  
   - Required fields: `sample_id`, `final_scale_intensity`  
   - Encoding: UTF‑8  

2. **Patient information file (.csv)**
   - Contains sample‑to‑patient group mapping  
   - Required fields: `Sample`, `patients_repeat`  
   - Encoding: UTF‑8  

3. **Glycan type file (.csv)**
   - Glycan composition and type mapping  
   - Required fields: `glycan`, `glycan.type`  
   - Encoding: UTF‑8  

### Files for Quantitative Analysis

1. **Missing‑value‑filled_zscore.csv**
   - Preprocessed expression matrix  

2. **sample_group_PCA.csv**
   - Sample grouping information  

## Analysis Module Details

### Glycopeptide Quantitative Submodules

1. **Heatmap Analysis**
   - CRLM vs CRC heatmaps  
   - CRLM vs LIHC heatmaps  
   - Combined gene‑screening heatmaps  
   - Z‑score normalized heatmaps  

2. **Regional Heatmap Analysis**
   - Mean‑ratio scatterplots  
   - Regional statistics  
   - Differential and z‑score heatmaps  

3. **Venn Diagram Analysis**
   - Comparison of different screening methods  
   - Overlap of significantly different genes  

4. **LDA & GO Analysis**
   - PCA, LDA, SVM  
   - GO enrichment  
   - Glycan type distribution  

5. **KEGG Analysis**
   - KEGG pathway enrichment  

6. **ROC Analysis**
   - Combined‑gene ROC curves  
   - Single‑gene ROC  

7. **Volcano Plot Analysis**
   - Differential expression visualization  

## Result Interpretation

### Common Visualizations

- **Heatmap**: Red = high expression, Blue = low expression  
- **Volcano plot**: X = log2 fold change, Y = ‑log10(p value)  
- **ROC curve**: Larger AUC = better diagnostic performance  
- **Venn diagram**: Shows overlap between gene sets  
- **GO/KEGG plots**: Show enriched biological pathways  

## FAQ

### Q1: What if sample matching fails?
Check sample IDs for exact consistency (case, spacing, etc.).

### Q2: Why is there "no enrichment result"?
Possible reasons:  
1) Too few significant genes  
2) Network issues  
3) Gene ID conversion failed  

### Q3: The analysis is slow?
Large datasets take time; close other programs and wait patiently.

### Q4: How to choose a module?
- Glycan analysis → glycan‑level research  
- Qualitative analysis → identification analysis  
- Quantitative analysis → differential expression research  

## Technical Support

- **Version**: v2.0  
- **Language**: R + Shiny  
- **Last updated**: 2024  

### Changelog (v2.0)

- Added glycopeptide quantitative module  
- Integrated 8 submodules  
- Enhanced visualization  
- Improved UI & stability  

### Usage Tips

- Read documentation before first use  
- For large datasets, consider batch processing  
- Download results promptly  
- Back up data regularly  

## License

MIT License (see LICENSE file)

## Contribution

1. Fork the project  
2. Create a feature branch  
3. Commit changes  
4. Submit a pull request  

## Contact

- Issue submission  
- Email: 1903180396@qq.com

---

**Note**: This system is designed for FFPE glycopeptide omics analysis. Ensure your data meets the required formats before use.

