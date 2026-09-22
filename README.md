# Enhanced Detection of Age-Related and Cognitive Declines Using Automated Hippocampal-To-Ventricle Ratio in Alzheimer's Patients

This repository contains the analysis code, derived measurements and
manuscript sources for the study above, published in *Human Brain Mapping*
(2025). We trained and evaluated three automated methods for segmenting the
hippocampus (HC) and the surrounding CSF-filled spaces (temporal horns of the
lateral ventricles) — Multi-Atlas Label Fusion (MALF), Nonlinear Patch-Based
segmentation (NLPB) and a Convolutional Neural Network (CNN) — and compared
them, together with FreeSurfer, on baseline T1w MRIs of 1641 participants from
the Alzheimer's Disease Neuroimaging Initiative (ADNI), using manual
segmentations of 80 cognitively healthy individuals as the gold standard. The
code here reproduces the statistical analyses, tables and figures of the paper.

## Looking for the HVR tool?

**If you want to compute hippocampal-to-ventricle ratios on your own data, do
not use this repository.** Use **hvr-cnn**, the maintained, containerised
version of the CNN segmentation method described in the paper:

- Source code, documentation and user guide:
  <https://github.com/soffiafdz/hvr-cnn>
- Container images:
  - `ghcr.io/soffiafdz/hvr-cnn`
  - `docker.io/soffiafdz/hvr_cnn` (the Docker Hub repository cited in the
    paper, updated to the new version)

hvr-cnn takes a T1w MRI (MINC or NIfTI — raw/native, already in stereotaxic
space, or AssemblyNet output) and produces hippocampus and temporal-horn
segmentations, their volumes, the hippocampal-to-ventricle ratio
(HVR = HC / (HC + VC)) and a QC image. It runs with Docker or podman, and with
Apptainer/Singularity on HPC systems. It is distributed under the GPL-3.0
licence.

The `container/` directory in *this* repository is the original, frozen
container source used for the paper. It is superseded by hvr-cnn and is kept
here only for the record; please use hvr-cnn instead.

## Citation

If you use this code or the HVR method, please cite:

> Fernandez-Lozano, S., Fonov, V., Schoemaker, D., Pruessner, J., Potvin, O.,
> Duchesne, S., Collins, D. L., & the Alzheimer's Disease Neuroimaging
> Initiative (2025). Enhanced detection of age-related and cognitive declines
> using automated hippocampal-to-ventricle ratio in Alzheimer's patients.
> *Human Brain Mapping*, 46(11), e70265. https://doi.org/10.1002/hbm.70265

```bibtex
@article{FernandezLozano2025HVR,
  author  = {Fernandez-Lozano, Sofia and Fonov, Vladimir and
             Schoemaker, Dorothee and Pruessner, Jens and
             Potvin, Olivier and Duchesne, Simon and
             Collins, D. Louis and
             {Alzheimer's Disease Neuroimaging Initiative}},
  title   = {Enhanced Detection of Age-Related and Cognitive Declines Using
             Automated Hippocampal-To-Ventricle Ratio in Alzheimer's
             Patients},
  journal = {Human Brain Mapping},
  year    = {2025},
  volume  = {46},
  number  = {11},
  pages   = {e70265},
  doi     = {10.1002/hbm.70265}
}
```

The article is open access (CC BY-NC-ND).

## Data availability

The MRI, demographic and clinical data used in this study come from the
Alzheimer's Disease Neuroimaging Initiative (ADNI) and **cannot be
redistributed here**. They are available from ADNI at
<https://adni.loni.usc.edu/>, subject to ADNI's application and data use
agreement.

What is and is not in this repository:

- **Included**: all code, the manuscript sources, and result objects that
  hold no participant-level information (effect sizes, correlation and
  permutation results under `data/rds/`).
- **Not included**: any image data; the ADNI tables that the parsing
  scripts read (`data/ADNIMERGE.csv`, `data/MRILIST.csv`,
  `data/UCSFFSX_11_02_15_20Nov2023.csv`,
  `data/UCSFFSX51_11_08_19_20Nov2023.csv`); and the **per-participant
  derived measurements** produced in this study (segmentation volumes,
  overlap values, QC decisions, scanner information, subject lists), which
  are keyed to ADNI participant identifiers and therefore fall under the
  ADNI Data Use Agreement. ADNI-approved investigators can obtain them from
  the corresponding author.
- **Not included**: the atlas/label libraries required by MALF and NLPB, and
  the CNN training set.

To re-run the analyses you therefore need an approved ADNI account, the
ADNI tables above placed in `data/`, and the derived measurements placed in
`data/derivatives/` and `lists/` as the scripts expect (file names are
given in the scripts).

## Repository layout

- `code/` — R, shell and Perl code (see
  [Contents of `code/`](#contents-of-code)).
- `data/` — tabular data (see [Data availability](#data-availability) for
  what is and is not included).
  - `data/derivatives/` — where the per-participant volumes and overlap
    metrics go (not included).
  - `data/rds/` — R objects produced by the analysis scripts; only those
    without participant-level rows are included.
- `lists/` — where subject lists, QC ratings and scanner information go
  (not included).
- `libraries/cnn/` — weights of the trained CNN ensembles
  (`ensemble_hcvc.pth`, `ensemble_hcvc-ag.pth`).
- `container/` — legacy container source for the CNN segmentation
  (superseded by [hvr-cnn](https://github.com/soffiafdz/hvr-cnn); see above).
- `quarto/` — Quarto sources of the manuscript (`manuscript.qmd`,
  `_quarto.yml`, `references.bib`, `styles.css`) and the final figures under
  `quarto/figures/`.
- `supplementaries/` — LaTeX source and compiled PDF of the supplementary
  material.
- `paper_submission/` — manuscript, cover letters, figures, tables and
  response to reviewers as submitted (`hbm1/`: first submission, `hbm2/`:
  revision).
- `environment.yml` — explicit conda specification (linux-64) of the
  environment used for the analyses.
- `renv.lock`, `renv/`, `.Rprofile` — R package environment.
- `LICENSE` — GPL-3.0.

Outputs written by the analysis scripts (`plots/`, `tables/`) are not tracked
in this repository; the final figures and tables as published are in
`quarto/figures/`, `supplementaries/` and `paper_submission/`.

## Reproducing the analyses

### Environments

Two environments are involved:

1. **conda** (`environment.yml`) — the full environment used for the
   published analyses. It is an explicit package specification for linux-64
   and includes MINC Toolkit v2 1.9.18.3, `minc2-simple`, PyTorch 1.10.1
   (CUDA 11.3), Python 3.9.13, Perl 5.32.1 and R 4.2.2. It is required only
   for the segmentation, QC-image and volume-extraction steps, which operate
   on MINC images.

   ```sh
   conda create --name hvr_validation --file environment.yml
   ```

2. **renv** (`renv.lock`) — the R packages used by the analysis scripts.
   `.Rprofile` activates renv automatically; from the repository root, start R
   and run:

   ```r
   renv::restore()
   ```

   The statistical analyses reported in the paper were run with **R 4.2.2**.
   The `renv.lock` currently in the repository was refreshed afterwards and
   records **R 4.5.2** with 168 packages; use `environment.yml` if you need
   the exact R version of the published run.

### Running the scripts

The R scripts locate files with the `here` package, so run them from the
repository root, for example:

```sh
Rscript code/data_parsing/parse_adnimerge-bl.R
```

Each script checks for the files it needs and sources its upstream script when
an intermediate object is missing, so running an analysis script will pull in
its dependencies. A full run follows the order below: data parsing and
cleaning, HVR calculation, then the analyses.

**Note**: the shell scripts under `code/backend/` and `code/segmentation/`
locate this repository automatically and take the locations of external
data from environment variables (`ADNI_PREPROC_DIR`, `MALF_HVR_DIR`,
`SNIPE_HVR_DIR`, `HVR_ADNI_DIR`, `MALF_LIBRARY_DIR`, `ICBM_HC_MODELS_DIR`);
each script says which one it needs.

## Contents of `code/`

### Data parsing and cleaning — `code/data_parsing/`

These scripts parse and clean the data used by the analyses and normally need
to be run first.

- `parse_adnimerge-bl.R` — parse the demographic data for the baseline
  images. Requires `data/ADNIMERGE.csv` and `data/MRILIST.csv` from ADNI
  (not included).
- `qc_segmentations_adni-bl.R` — filter out images and segmentations that
  failed QC. Needs the QC lists and segmentation volumes (not included, see
  above).
- `parse_freesurfer-vols.R` — parse ADNI's FreeSurfer results and our own
  FreeSurfer run. Requires `data/UCSFFSX_11_02_15_20Nov2023.csv` and
  `data/UCSFFSX51_11_08_19_20Nov2023.csv` from ADNI (not included).
- `extract_demog_adni-bl.R` — build the demographics table. Requires the
  scripts above to have been run.

### HVR calculation and analyses — `code/analysis/`

- `adjust_hc-hvr_adni-bl.R` — adjust for head size and compute HVR for all
  segmentation methods. Requires `data/ADNIMERGE.csv` from ADNI (not
  included).
- `compare_man-seg.R` — compare the segmentation methods against the manual
  labels: overlap similarity, correlation of volumetry and Bland-Altman
  plots.
- `compare_adni-bl.R` — compare the segmentation methods on the ADNI data:
  table of failures, table of HC volume and HVR, effect sizes between
  cognitively healthy participants and patients with AD, and cross-correlation
  of volumetry between methods.
- `correlate_age-memory_adni-bl.R` — associations of HC volume and HVR with
  age, RAVLT learning scores and ADAS13.

### Automatic segmentation — `code/segmentation/`

#### Multi-Atlas Label Fusion (MALF)

Source code: `code/segmentation/source_code/malf/` (`seg_hippo.pl`,
`local_hc_reg`, `pp_hc_xcorr_only`, `minctracc-w-nmi`).

**NOTE**: a library of labels and T1w images with their transformations to a
standard template is required to run this algorithm and is not included in
this repository.

Application:

- Manually segmented images:
  `code/segmentation/application/manual-seg_malf_hcvc.sh`
- ADNI (baseline): `code/segmentation/application/adni-bl_malf_hcvc.sh`

#### Nonlinear Patch-Based segmentation (NLPB)

Source code: `code/segmentation/source_code/nlpb/snipe_minipipe_reduc.pl`.

**NOTE**: a library of labels and T1w images, as well as the ROI masks, are
required to run this algorithm and are not included in this repository.

Application:

- Manually segmented images:
  `code/segmentation/application/manual-seg_nlpb_hcvc.sh`
- ADNI (baseline): `code/segmentation/application/adni-bl_nlpb_hcvc.sh`

#### Convolutional Neural Network (CNN)

The network was trained with the `py_deep_seg` framework (V. S. Fonov),
which is not included in this repository. The
inference code used for the paper is in `container/`, and the maintained
version of it is [hvr-cnn](https://github.com/soffiafdz/hvr-cnn).

**NOTE**: the training library is not included in this repository. The final
weights of the trained networks (5-fold cross-validation) are in
`libraries/cnn/ensemble_hcvc.pth` (HC and CSF) and
`libraries/cnn/ensemble_hcvc-ag.pth` (HC head/body/tail, CSF and amygdala).

Application:

- Manually segmented images:
  - `code/segmentation/application/manual-seg_cnn_hcvc_hcvc-ag.sh`
  - `code/segmentation/application/manual-seg_cnn_hcvc_validation.sh`
- ADNI (baseline):
  `code/segmentation/application/adni-bl_cnn_hcvc_hcvc-ag.sh`

### QC image generation — `code/qc/`

QC images were created to visually evaluate the segmentations.

- `qc_plot.pl` — QC images of the segmentation labels (separate head, body
  and tail regions).
- `qc_plot_reduc.pl` — QC images of the segmentation labels (whole
  hippocampus).

### Backend scripts — `code/backend/`

These shell and R scripts supported the analysis (data management, volume
extraction, document assembly) but are not required to reproduce the results.

- `get_brainmask_qc-images.sh` — obtain brain mask QC images.
- `compare_man-seg.sh` — overlap similarity between the MALF, NLPB and CNN
  segmentations and the manual labels.
- `compare_man-seg_validation.sh` — the same comparison on the validation
  datasets (ADNI and ICBM).
- `calculate_kappa_man-seg_fs.sh` — overlap similarity between the FreeSurfer
  segmentations and the manual labels.
- `compare_cnns_adni-bl.sh` — overlap similarity between the two trained CNNs
  (whole HC versus head/body/tail).
- `create_list_adni_baseline.sh` — extract the first recorded session of each
  ADNI participant.
- `create_list_adni_preproc.sh` — extract all preprocessed sessions of each
  ADNI participant.
- `extract_icc_scale-factor_adni-bl.sh` — extract the intracranial cavity
  volume (ICC) and the stereotaxic scale factor from the preprocessed ADNI
  data.
- `extract_scanner_adni-bl.sh` — extract scanner field strength from the MRI
  headers.
- `extract_volumes_hcvc.sh` — compute volumes from the segmentations
  (HC and CSF).
- `extract_volumes_hcvc-ag.sh` — compute volumes from the segmentations
  (HC and CSF head/body/tail, and amygdala).
- `fill_qc-dirs_adni-bl.sh` — create directories of QC images for curation.
- `link_files_adni-bl.sh` — create symlinks of the data.
- `relabel_cnn_hcvc-ag_hcvc.sh` — relabel head/body/tail segmentations to
  whole HC and CSF.
- `resample_segmentations_adni-bl_nlpb.sh` — resample the NLPB segmentations
  back to the standard space.
- `compile_figures.sh`, `compile_figures_txt.sh`, `compile_figures_word.R` —
  assemble the figures into a single LaTeX or Word document.
- `compile_tables.sh`, `compile_tables_word.R` — assemble the tables into a
  single LaTeX or Word document.
- `compile_supplementaries.sh` — assemble and compile the supplementary
  material.

## Licence

This repository is released under the GNU General Public License v3.0; see
`LICENSE`.

The legacy `container/model/` directory contains only the network code
needed for inference (by V. S. Fonov, included with permission).

## Contact

Questions about the code are best raised as issues on this repository
(for the HVR tool itself, use the
[hvr-cnn](https://github.com/soffiafdz/hvr-cnn) issue tracker).

Corresponding author: D. Louis Collins (louis.collins@mcgill.ca), McConnell
Brain Imaging Centre, Montreal Neurological Institute, McGill University.

## Acknowledgements

This project was made possible by the Brain Canada Foundation, through the
Canada Brain Research Fund, with the financial support of Health Canada, the
Canadian Institutes of Health Research Project Grant FRN 165921, and
La Fondation Famille Louise & André Charron.

Data used in the preparation of this work were obtained from the Alzheimer's
Disease Neuroimaging Initiative (ADNI) database (<https://adni.loni.usc.edu/>).
As such, the investigators within the ADNI contributed to the design and
implementation of ADNI and/or provided data but did not participate in the
analysis or writing of this report. A complete listing of ADNI investigators
is available at
<http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf>.
