# Data Format Examples

This document provides data format examples for the input files required by the glycopeptide analysis system.

## 1. Inter-group Normalization Data File (zujian_normalized_data.csv)

### Required Columns:
- `Glycopeptides`: Glycopeptide sequence
- `Protein.Accession`: Protein accession number
- `raw.file`: Raw file name
- `sample`: Sample name
- `final_scale_intensity`: Final scaled intensity value (log2 transformed value)
- `sample_id`: Sample ID

### Example Data:
```csv
Glycopeptides,Protein.Accession,raw.file,sample,final_scale_intensity,sample_id
NLTEEQIAEFK-N2H5F1S1G0,P12345,FFPE_mix1_glyco_4h_3ug_60K_repeat1.raw,FFPE_mix1_glyco_4h_3ug_60K_repeat1,15.234,FFPE_mix1_glyco_4h_3ug_60K_repeat1
NLTEEQIAEFK-N2H5F1S1G0,P12345,FFPE_mix1_glyco_4h_3ug_60K_repeat2.raw,FFPE_mix1_glyco_4h_3ug_60K_repeat2,14.876,FFPE_mix1_glyco_4h_3ug_60K_repeat2
NLTEEQIAEFK-N2H5F1S1G0,P12345,FFPE_mix1_glyco_4h_3ug_60K_repeat3.raw,FFPE_mix1_glyco_4h_3ug_60K_repeat3,15.123,FFPE_mix1_glyco_4h_3ug_60K_repeat3
VYVEELKPTPEGDLEILLQK-N2H5F1S1G0,P67890,FFPE_mix1_glyco_4h_3ug_60K_repeat1.raw,FFPE_mix1_glyco_4h_3ug_60K_repeat1,12.456,FFPE_mix1_glyco_4h_3ug_60K_repeat1
VYVEELKPTPEGDLEILLQK-N2H5F1S1G0,P67890,FFPE_mix1_glyco_4h_3ug_60K_repeat2.raw,FFPE_mix1_glyco_4h_3ug_60K_repeat2,12.789,FFPE_mix1_glyco_4h_3ug_60K_repeat2
```

### Notes:
- `final_scale_intensity` should be log2 transformed values
- `sample_id` must exactly match the `Sample` column in the patient information file
- Glycopeptide sequence format: `peptide_sequence-glycan_composition`

## 2. Patient Information File (FFPE_patients_reflect.csv)

### Required Columns:
- `Sample`: Sample ID (corresponds to sample_id in inter-group normalization data)
- `patients_repeat`: Patient repeat information

### Example Data:
```csv
Sample,patients_repeat
FFPE_mix1_glyco_4h_3ug_60K_repeat1,CRC.1_repeat1
FFPE_mix1_glyco_4h_3ug_60K_repeat2,CRC.1_repeat2
FFPE_mix1_glyco_4h_3ug_60K_repeat3,CRC.1_repeat3
FFPE_mix2_glyco_4h_3ug_60K_repeat1,CRC.2_repeat1
FFPE_mix2_glyco_4h_3ug_60K_repeat2,CRC.2_repeat2
FFPE_mix2_glyco_4h_3ug_60K_repeat3,CRC.2_repeat3
FFPE_mix3_glyco_4h_3ug_60K_repeat1,CRLM.1_repeat1
FFPE_mix3_glyco_4h_3ug_60K_repeat2,CRLM.1_repeat2
FFPE_mix3_glyco_4h_3ug_60K_repeat3,CRLM.1_repeat3
FFPE_mix4_glyco_4h_3ug_60K_repeat1,LIHC.1_repeat1
FFPE_mix4_glyco_4h_3ug_60K_repeat2,LIHC.1_repeat2
FFPE_mix4_glyco_4h_3ug_60K_repeat3,LIHC.1_repeat3
```

### Format Description:
- `patients_repeat` format: `patient_name.cancer_type_repeat_number`
- Supported cancer types: CRC (Colorectal Cancer), CRLM (Colorectal Cancer Liver Metastasis), LIHC (Liver Hepatocellular Carcinoma)
- The system will automatically parse cancer type and repeat number

## 3. Glycan Type File (glycan-glycantype.csv)

### Required Columns:
- `glycan`: Glycan sequence
- `glycan.type`: Glycan type

### Example Data:
```csv
glycan,glycan.type
N2H5F1S1G0,Complex
N2H5F1S0G0,Complex
N2H5F0S1G0,Complex
N2H5F0S0G0,Complex
N2H4F1S1G0,Complex
N2H4F1S0G0,Complex
N2H4F0S1G0,Complex
N2H4F0S0G0,Complex
N2H3F1S1G0,Hybrid
N2H3F1S0G0,Hybrid
N2H3F0S1G0,Hybrid
N2H3F0S0G0,Hybrid
N2H2F1S1G0,High-mannose
N2H2F1S0G0,High-mannose
N2H2F0S1G0,High-mannose
N2H2F0S0G0,High-mannose
```

### Glycan Composition Description:
- N: N-acetylglucosamine (GlcNAc)
- H: Mannose
- F: Fucose
- S: Sialic Acid
- G: Glucose

### Glycan Types:
- Complex: Complex-type glycans
- Hybrid: Hybrid-type glycans
- High-mannose: High-mannose type glycans

## Data Preparation Recommendations

### 1. Data Quality Control
- Check missing value ratio
- Ensure sample ID matching
- Verify data format correctness

### 2. File Encoding
- Recommend using UTF-8 encoding
- Avoid using special characters

### 3. Data Scale
- Supports large files (maximum 1000MB)
- Recommend batch processing for very large datasets

### 4. Backup Recommendations
- Keep original data backup
- Record data processing steps

## Frequently Asked Questions

### Q: Why does an error appear after uploading files?
A: Please check:
- Whether the file format is CSV
- Whether column names match exactly
- Whether file encoding is UTF-8
- Whether sample IDs match

### Q: How to handle missing values?
A: The system will automatically use KNN algorithm to impute missing values, but recommend:
- Check the cause of missing values
- Ensure missing value ratio does not exceed 50%
- Consider data preprocessing steps

### Q: What cancer types are supported?
A: Currently supported:
- CRC (Colorectal Cancer)
- CRLM (Colorectal Cancer Liver Metastasis)
- LIHC (Liver Hepatocellular Carcinoma)

If you need to add other types, please modify the relevant parts in the code.