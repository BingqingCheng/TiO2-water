# Data repository for "Mechanistic insight on water dissociation on pristine low-index TiO$_2$ surfaces from machine learning molecular dynamics simulations"

arXiv preprint: \url{https://arxiv.org/abs/2303.07433}

## MLPs:
Contains machine learning potentials trained on optb88vdw (./optb88vdw/), and the Delta-learning potentials (./delta-pbe/, ./delta-scan-1200Ry/, ./delta-scan-350Ry/). The Delta-learning potentials need to be used together with the optb88vdw MLP.

## cp2k-input:
The DFT input files for CP2K.

## fig1-data:
Source data for generating Fig.1. Contains the free energy surfaces, water density profiles, and water orientations.

## metad-input:
The example input files for running metadynamics simulations using LAMMPS and PLUMED.

## training-data:
The training data for the optb88vdw MLPs, in N2P2 format. Units are in hartree/Bohr.

## delta-PBE-training-data:
The training data for the Delta-learning PBE  MLPs, in N2P2 format. Units are in hartree/Bohr.

## delta-SCAN-training-data:
The training data for the Delta-learning SCAN MLPs, in N2P2 format. Units are in hartree/Bohr.
