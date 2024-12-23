# README: Generating Chord Plots for Brain-EMG Connectivity

This README provides step-by-step instructions on how the R script generates chord plots to visualize brain-EMG connectivity data. The process involves reading connectivity matrices, performing statistical tests, and creating plots with the `circlize` package.

## Overview

### Key Steps:
1. Load connectivity data from `.mat` files.
2. Perform statistical tests to identify significant connections.
3. Mask non-significant connections.
4. Rearrange connectivity data for proper visualization.
5. Generate and customize chord plots.

## Step-by-Step Process

### 1. Load Connectivity Data

- **Input Files**: Connectivity matrices (`netVals_young_...`) stored as `.mat` files.
- **Tools Used**: `R.matlab` library for reading MATLAB `.mat` files.
- **Purpose**:
  - Load subject-wise and average connectivity data for specific trial conditions (e.g., `"LEI"`, `"LME"`) and frequency bands (e.g., `"theta"`, `"alpha"`).

### 2. Perform Statistical Tests

- **Goal**: Identify significant connections between components.
- **Approach**:
  - Perform t-tests on connectivity values using `t.test()` for each pair of components (excluding self-connections).
  - Reshape the results into matrices (`pvals_theta_box` and `pvals_alpha_box`).
  - Apply corrections (e.g., FDR) if necessary.

### 3. Mask Non-Significant Connections

- **Method**:
  - Use significance thresholds (`p < 0.05`) to create binary masks (`isSig_theta_box`, `isSig_alpha_box`).
  - Apply masks to average connectivity data (`thetaAve`, `alphaAve`) to zero out non-significant connections.

### 4. Rearrange Data for Visualization

- **Steps**:
  - Rearrange connectivity matrices to a specific order for better visual interpretation.
  - Separate data into cortico-cortical and cortico-muscular connections.

### 5. Generate Chord Plots

- **Library**: `circlize`
- **Customizations**:
  - Define component labels (`factors`) and colors (`color_vals`).
  - Adjust plot parameters (`circos.par`) to optimize the layout.
  - Add links to represent significant connections:
    - Red for positive correlations.
    - Blue for negative correlations.
    - Gray for non-significant connections.

- **Output**: Save plots as `.png` files.

## Example Workflow

### Initialization
```r
library(R.matlab)
library(circlize)
```

### Read Data
```r
netData_sbjs <- readMat("netVals_young_LEI_early_epoch_PPT_aveTime_sbjs.mat")
```

### Statistical testing
```r
pvals_theta <- t.test(netData_sbjs$netVals[[1]][[1]])$p.value
pvals_alpha <- t.test(netData_sbjs$netVals[[2]][[1]])$p.value
```

### Mask and Rearrange Data
```r
thetaAve <- netData_ave$netVals.ave[[1]]
alphaAve <- netData_ave$netVals.ave[[2]]
thetaAve[isSig_theta_box == 0] <- 0
alphaAve[isSig_alpha_box == 0] <- 0
```

### Create Chord Plot
```r
circos.initialize(factors, xlim = cbind(c(0, 0), 2))
circos.trackPlotRegion(ylim = c(0, 1))
circos.link("ACC", c(0, 0.5), "SMA_r", c(0, 0.3), col = "red")
```
```mermaid
flowchart TD
    A[Load Connectivity Data] --> B[Perform Statistical Tests]
    B --> C[Mask Non-Significant Connections]
    C --> D[Rearrange Data for Visualization]
    D --> E[Generate Chord Plots]
    E --> F[Save Plots as PNG Files]
