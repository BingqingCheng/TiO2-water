# Data repository for "Mechanistic insight on water dissociation on pristine low-index TiO$_2$ surfaces from machine learning molecular dynamics simulations"

[![DOI](https://zenodo.org/badge/613786300.svg)](https://zenodo.org/badge/latestdoi/613786300)

arXiv preprint: https://arxiv.org/abs/2303.07433

## MLPs:
Contains machine learning potentials trained on optb88vdw (./optb88vdw/), and the Delta-learning potentials (./delta-pbe/, ./delta-scan-1200Ry/, ./delta-scan-350Ry/). The Delta-learning potentials need to be used together with the optb88vdw MLP.

## cp2k-input:
The DFT input files for CP2K.

## fig1-data:
Source data for generating Fig.1. Contains the free energy surfaces, water density profiles, and water orientations.

## fig2-data:
Source data for generating the kPCA maps in Fig.2.

## fig3-data:
Source data for generating the kPCA maps in Fig.3a,b, and the transition matrix in Fig.3d.

## metad-input:
The example input files for running metadynamics simulations using LAMMPS and PLUMED.

## training-data:
The training data for the optb88vdw MLPs, in N2P2 format. Units are in hartree/Bohr.

## delta-PBE-training-data:
The training data for the Delta-learning PBE  MLPs, in N2P2 format. Units are in hartree/Bohr.

## delta-SCAN-training-data:
The training data for the Delta-learning SCAN MLPs, in N2P2 format. Units are in hartree/Bohr.

## example-analysis:
Python script and notebook to generate features for hydrogen environments, perform classification and kernal PCA.
