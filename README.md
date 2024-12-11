# Step 6: Statistical Significance Tests for Corticomuscular Connectivity

This step focuses on performing statistical tests to identify significant brain-to-muscle and muscle-to-brain connectivity from the connectivity matrices derived in earlier steps. The analysis includes comparing connectivity values against significance thresholds and visualizing the results.

## Table of Contents
1. [Initial Setup](#initial-setup)
2. [Parameters and Configuration](#parameters-and-configuration)
3. [Subject and Group Selection](#subject-and-group-selection)
4. [Defining Phases and Trials](#defining-phases-and-trials)
5. [Running the Statistical Analysis](#running-the-statistical-analysis)
6. [Output and Visualization](#output-and-visualization)

---

## 1. Initial Setup

```matlab
clearvars -except subj groupName earlyFraction;
close all; clc;
```

- **`clearvars`**: Clears all variables except for `subj`, `groupName`, and `earlyFraction` if they exist.
- **`close all`**: Closes all figure windows.
- **`clc`**: Clears the command window.

---

## 2. Parameters and Configuration

```matlab
fs = string(filesep) + string(filesep);
fPath = string(pwd) + fs;
global timewindow;
timewindow = [-0.4, -0.1];
alpha = 0.01; % Significance level for statistical tests
windowLengthSec = 0.4;
windowStepSizeSec = 0.02;
GUI_MODE = 'nogui';
VERBOSITY_LEVEL = 0;
```

- **`timewindow`**: Time window for analysis (e.g., from -400 ms to -100 ms).
- **`alpha`**: Significance threshold (0.01 in this case).
- **`windowLengthSec`**: Length of the sliding window in seconds.
- **`windowStepSizeSec`**: Step size for the sliding window.
- **`GUI_MODE`**: Whether to show GUI (`'nogui'` disables it).
- **`VERBOSITY_LEVEL`**: Controls the level of output detail.

---

## 3. Subject and Group Selection

```matlab
if ~exist('groupName', 'var') || isempty(groupName)
    groupName = "young";
else
    groupName = string(groupName);
end

if ~exist('earlyFraction', 'var') || isempty(earlyFraction)
    earlyFraction = 0.999;
end

if groupName == "young"
    subjs = ["PS04", "PS05", ..., "PS27"];
elseif groupName == "old"
    subjs = ["PS08", "PS09", ..., "PS36"];
end
```

- **`groupName`**: Specifies the group (`"young"` or `"old"`).
- **`earlyFraction`**: Fraction of perturbations to categorize as early.

---

## 4. Defining Phases and Trials

```matlab
trialTypes = ["LEI", "LME"];
phases = ["early_epoch_PPT"];
```

- **`trialTypes`**: Types of trials (e.g., `"LEI"` for left extension initiation).
- **`phases`**: Defines which phases to analyze based on the `earlyFraction`.

---

## 5. Running the Statistical Analysis

The script runs statistical tests to identify significant connections based on predefined thresholds and conditions.

- **Connectivity Matrix Processing**: Processes the matrices to compute brain-to-muscle and muscle-to-brain connectivity.
- **Statistical Comparisons**: Tests connectivity against surrogate or null distributions to identify significant connections.

---

## 6. Output and Visualization

- **Significant Connections**: Outputs significant brain-to-muscle and muscle-to-brain connections.
- **Figures**: Generates plots to visualize significant connectivity patterns.

---

## Notes

- Ensure all dependencies (EEGLAB, DIPFIT, etc.) are installed.
- Modify `alpha` for different significance thresholds.
- Adjust `windowLengthSec` and `windowStepSizeSec` for different temporal resolutions.
