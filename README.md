# Spatiotemporal-multimodal integration reveals a BCG-induced skin-blood crosstalk

*reproduce all analyses using this codebase*

---

## Table of Contents

- [Overview](#overview)  
- [Getting Started](#getting-started)  
- [Running the Scripts](#running-the-scripts)  
- [How to Cite](#how-to-cite)  

---

## Overview

This repository contains scripts for data exploration, differential expression analysis, and multi-modal data integration. The project uses R and R Markdown to process, analyze, and visualize results from diverse data sources, including cfRNA, dermatological images, and GeoMx spatial transcriptomics data.

---

## Getting Started

### 1) clone repo

```
git clone https://github.com/CompBio-Lab/bcg_skinblood_crosstalk.git
```

### 2) get data
> download data from this [zenodo link]()
> put data in bcg_skinblood_crosstalk folder

### 3) open .Rproj file
> bcg_skinblood_crosstalk.Rproj

### 4) Install renv if you don't have it:

```
install.packages("renv")
```

### 5) Restore the project environment from renv.lock:

```
renv::restore()
```
> This will install the exact package versions used when this project was last saved.


## Running scripts

- Open .Rmd files in RStudio

- Knit the files to HTML or run interactively

- Explore the folder structure above to locate analyses of interest


## How to Cite

- [preprint]()