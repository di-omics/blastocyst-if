# blastocyst-if

3D immunofluorescence analysis pipeline for whole-mount preimplantation blastocysts.

## Stack
- Cellpose / Cellpose-SAM (3D nuclear segmentation, MPS-accelerated on Apple Silicon)
- napari (interactive 3D review)
- bioio (CZI / LIF / OME-TIFF readers)
- scikit-image, pandas (measurements, export)

## Setup
```bash
conda env create -f environment.yml -p /opt/shared-envs/blastocyst-if
conda activate /opt/shared-envs/blastocyst-if
```

## Starting dataset
Niakan lab, Simon et al. 2025 (Nat Commun, ERK suppression in human blastocyst).
Confocal stacks: figshare 28597145.
Reference pipeline: zenodo 15640446.
