# Frobenius–Schur Inversion: Numerical Experiments in MATLAB

This repository contains the MATLAB experiments conducted for the master's thesis: 
**"Numerical methods for inverting complex matrices: A focus on Frobenius-Schur inversion."**

## Overview
The experiments evaluate the **accuracy** and **efficiency** of different numerical methods for computing the inverse of complex matrices.  
In particular, the repository mainly compares:
- **Frobenius-Schur inversion method**  
- **Default LU-decomposition-based inversion method** (MATLAB built-in)
- **QR-decomposition-based inversion method and SVD-decomposition-based inversion method**

## Repository Structure
- '.m' files in the **root folder**:  
  Each script named after 'Test_xxx' corresponds to one experiment. Running a script will generate the associated '.png' figure(s) that illustrate the results.  
- '.png' files:  
  Precomputed results obtained by running the corresponding '.m' script (same filename).  
- 'help_functions/':  
  Contains supporting MATLAB functions used across experiments.

## How to Reproduce the Experiments
1. Open MATLAB.  
2. Navigate to the project root folder.  
3. Run any of the main '.m' files (e.g., 'Test_5_1_1_1.m', 'Test_5_1_1_2.m', ...).  
4. The script will generate the corresponding .png figures, reproducing the results shown in the thesis.

## Requirements
MATLAB (R2022a or later)

