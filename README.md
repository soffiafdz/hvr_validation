# Automatization and validation of the hippocampal-to-ventricle ratio in a clinical sample

This repository contains the code to reproduce the analyses, tables and figures
for the manuscript.

Source MRI and demographic data can be requested to ADNI.

## Contents:
- `lists`: Directory with data in the form of lists and CSV files.
- `data/`: Directory with output data.
    - `data/derivatives:` Subirectory for derivatives created from the analyses.
    - `data/rds`: Subdirectory with RDS objects created from the analyses.
- `libraries/cnn`: Weigths for the trained CNN.
- `code/`: Directory with needed R and shell code.
    - `code/segmentation`: Subdirectory with the code used to run the automatic
    segmentation algorithms.
    - `code/qc`: Scripts used to generate the QC images that were used to
    evaluate the output of the automatic segmentations.
    - `code/data_parsing`: Code to clean and parse the data needed for the
    analyses.
    - `code/analysis`: Code needed for the analyses.
    - `code/backend`: Subdirectory with other code used for managing data
    but unnecessary to reproduce the results.
- `renv/`: Renv sourcedirectory for replicating R library.
- renv.lock: Renv lock file.
- LICENSE: GPL-3.0 license.
- README.md: This file.

## Running automatic segmentation algorithms

### Multi-atlas Label Fusion

#### Source code

The source code to run MALF can be found on the
`code/segmentation/source_code/malf/` directory.

**NOTE**: A library of labels and T1w images and their respective
transformation files to a standard template is needed to run this algorithm.
This library is not included in this repository.

#### Application

The following scripts run the MALF algorithm on the data.
- Manually segmented images:
    - `code/segmentation/application/manual-seg_malf_hcvc.sh`
- ADNI (baseline):
    - `code/segmentation/application/adni-bl_malf_hcvc.sh`

### Non-local patch based segmentation

#### Source code

The pipeline to run NLPB is located on:
`code/segmentation/source_code/nlpb/snipe_minipipe_reduc.pl`.

**NOTE**: A library of labels and T1w images, as well as the masks for the
ROI are needed to run this algorithm.
This library is not included in this repository.

#### Application

The following scripts run NLPB on the data:
- Manually segmented images:
    - `code/segmentation/application/manual-seg_nlpb_hcvc.sh`
- ADNI (baseline):
    - `code/segmentation/application/adni-bl_nlpb_hcvc.sh`

### Convolutional Neural Network

#### Source code

The code for the CNN training and inference of HC and Amygdala segmentation is
on the `code/segmentation/source_code/py_deep_seg/` directory.

**NOTE**: The training library used for training the network is not included
in this repository. On the other hand, the final weights of the trained
network (5-CV) can be found on
`libraries/cnn/ensemble_hcvc.pth` and `libraries/cnn/ensemble_hcvc-ag.pth`.

#### Application

The following scripts run the CNN on the data:
- Manually segmented images:
    - `code/segmentation/application/manual-seg_cnn_hcvc_hcvc-ag.sh`
- ADNI (baseline):
    - `code/segmentation/application/adni-bl_cnn_hcvc_hcvc-ag.sh`

## QC image generation

QC images were created to evaluate the segmentations using these two scripts.

- `code/qc/qc_plot.pl`
    - Script for creating QC images of the segmentation labels
      (separate regions for head/body/tail).
- `code/qc/qc_plot_reduc.pl`
    - Script for creating QC images of the segmentation labels
      (whole hippocampus).

## Data parsing and cleaning

These scripts parse and/or clean the data to be used for the analysis.
They usually need to be run before the following sections.

- `code/data_parsing/parse_adnimerge-bl.R`
    - Parse needed demographic data from ADNI for the baseline images.
    - **NOTE**: the files `data/ADNIMERGE.csv` and `data/MRILIST.csv`
    are required to run this script but are not included in this repo.
    They can be acquired directly from ADNI with the appropriate permissions.
- `code/data_parsing/qc_segmentations_adni-bl.R`
    - Filter out images and/or segmentations that failed QC.
    - The QC lists and HC volumes from the segmentations required to run this
    script are all included in this repo.
- `code/data_parsing/parse_freesurfer-vols.R`
    - Parse ADNI's and our own implementations of FreeSurfer.
    - The required files:
    `data/UCSFFSX_11_02_15_20Nov2023.csv` and
    `data/UCSFFSX51_11_08_19_20Nov2023.csv`
    are required to run this script bur are not included in this repo.
    They can be acquired directly from ADNI with the appropriate permissions.
    - For reproducibility purposes, the output generated from this script is
    included in this repo.

## Data processing (HVR calculation)
- `code/analysis/adjust_hc-hvr_adni-bl.R`
    - Adjust for head size and calculate HVR for all segmentation methods.
    - **NOTE**: the file `data/ADNIMERGE.csv` is required to run this script
    but is not included in this repo.
    It can be acquired directly from ADNI with the appropriate permissions.
    - For reproducibility purposes, the output generated from this script is
    included in this repo.

## Demographics table
- `code/data_parsing/extract_demog_adni-bl.R`
    - Clean and create a table for demographic data
    - The previous sections need to be run first.

## Segmentation method comparisons
- `code/analysis/compare_cnns_adni-bl.R`
    - Compare the overlap similarity between two trained CNNs
    (whole HC and head/body/tail HC):
      - Plot (not included in the paper).
- `code/analysis/compare_freesurfer_adni-bl.R`
    - Compare the overlap similarity between the different versions of FreeSurfer:
      - Plot (not included in the paper).
- `code/analysis/compare_segmentations_manual-labels.R`
    - Compare the performance of the segmentation models:
        - Overlap similarity with manual labels.
        - Correlation of volumetry with manual labels.
        - Bland-altman plots.
- `code/analysis/compare_segmentations_adni-bl.R`
    - Compare the segmentation models on the ADNI data:
        - Table of failures.
        - Table of HCv and HVR calculations.
        - CH:AD effect sizes.
        - Similarity (cross-correlation) of volumetry between methods.

## HC associations with age, memory and global cognition
- `code/analysis/correlate_age-memory_adni-bl.R`
    - Evaluate the association of the calculated HCv and HVR with age, RAVLT and ADAS13.

## Backend scripts

These shell scripts were written to facilitate the analysis but are not required
to reproduce the results.

- `code/backend/get_brainmask_qc-images.sh`
    - Obtain Brain mask QC images.
- `code/backend/calculate_kappa_man-seg_cnn.sh`
    - Compare overlap similarity between CNN segmentations and manual labels.
- `code/backend/calculate_kappa_man-seg_fs.sh`
    - Compare overlap similarity between FreeSurfer segmentations and manual labels.
- `code/backend/compare_cnns_adni-bl.sh`
    - Compare the overlap similarity between two trained CNNs
    (whole HC and head/body/tail HC).
- `code/backend/create_list_adni_baseline.sh`
    - Extract the first recorded session of all ADNI participants.
- `code/backend/create_list_adni_preproc.sh`
    - Extract all the sessions that have been preprocessed of all ADNI participants.
- `code/backend/extract_icc_scale-factor_adni-bl.sh`
    - Extract ICC and Scale factor from the preprocessed ADNI participants.
- `code/backend/extract_scanner_adni-bl.sh`
    - Extract scanner fieldstrength data from the header of the MRI data.
- `code/backend/extract_volumes_hcvc-ag.sh`
    - Calculate volume from the segmentations (HC & CSF head/body/tail & Amygdala).
- `code/backend/extract_volumes_hcvc.sh`
    - Calculate volume from the segmentations (HC & CSF).
- `code/backend/fill_qc-dirs_adni-bl.sh`
    - Create directories of QC images for curation.
- `code/backend/link_files_adni-bl.sh`
    - Create symlinks of the data.
- `code/backend/relabel_cnn_hcvc-ag_hcvc.sh`
    - Relabel head/body/tail segmentations to full HC & CSF.
- `code/backend/resample_segmentations_adni-bl_nlpb.sh`
    - Resample NLPB segmentations to original standard space.
