# Let’s Explore and Learn to Analyze Untargeted Metabolomics Data

[![License: CC BY-NC
4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![DOI](https://zenodo.org/badge/849331979.svg)](https://doi.org/10.5281/zenodo.15062929)

## Welcome to **Metabonaut**! 🚀

Metabonaut presents a series of workflows based on a small LC-MS/MS
dataset, utilizing **R** and **Bioconductor** packages. These workflows
demonstrate how to adapt various algorithms to specific datasets and
seamlessly integrate R packages for efficient, reproducible data
processing.

## Available Vignettes

### 1. [Complete End-to-End LC-MS/MS Metabolomic Data Analysis](https://rformassspectrometry.github.io/Metabonaut/articles/a-end-to-end-untargeted-metabolomics.html)

This primary workflow guides you through each step of the analysis, from
**preprocessing raw data** to **statistical analysis** and **metabolite
annotation**. 📄 [View Full R Code
(Quarto)](https://github.com/rformassspectrometry/Metabonaut/blob/main/vignettes/a-end-to-end-untargeted-metabolomics.qmd)

### 2. [Dataset Investigation](https://rformassspectrometry.github.io/Metabonaut/articles/dataset-investigation.html)

Before diving into the analysis, learn about **key aspects to examine in
your dataset** to ensure smooth processing and prevent troubleshooting
issues later in the pipeline.

### 3. [Seamless Alignment: Merging New Data with an Existing Preprocessed Dataset](https://rformassspectrometry.github.io/Metabonaut/articles/alignment-to-external-dataset.html)

Discover how to use a **flexible alignment algorithm** to integrate new
datasets with previously processed ones based on features of interest.

### 4. [Quality Control and Feature Selection Using *notame*](https://rformassspectrometry.github.io/Metabonaut/articles/notame_normalization_and_feature_selection.html)

Explore **notame** as a robust alternative for normalization and feature
selection. This vignette covers essential steps for **quality control,
normalization, and feature selection** to ensure your data is clean and
ready for analysis.

### 5. [LC-MS/MS Data Annotation Using R and Python](https://rformassspectrometry.github.io/Metabonaut/articles/SpectriPy_tutorial_metabonaut.html)

Explore the **SpectriPy** package for **LC-MS/MS data annotation**. This
tutorial demonstrates how to combine the strengths of Python and R MS
libraries for comprehensive annotation.

### 6. [Large Scale Processing Using *xcms*](https://rformassspectrometry.github.io/Metabonaut/articles/large-scale-analysis.html)

While *xcms* is known for its scalability, this guide shows you how to
practically handle **large-scale dataset processing (\>4,000 files)** on
standard hardware.

### 7. [Using and Creating Metabolomics Data Annotation Resources](https://rformassspectrometry.github.io/Metabonaut/articles/creating-using-annotation-resources.html)

Learn how to **create and utilize custom metabolomics data annotation
resources** in R. This vignette guides you through building your own
annotation databases and integrating them into your analysis workflow.
We demonstrate by downloading data from GNPS and MassBank libraries.

------------------------------------------------------------------------

For a full list of all available vignettes, visit the [**Metabonaut
website**](https://rformassspectrometry.github.io/Metabonaut/).

------------------------------------------------------------------------

## 📌 Reproducibility & Updates

We strive for **reproducibility**. These workflows are designed to
remain stable over time, allowing you to run all vignettes together as
one comprehensive **super-vignette**.

- **Major updates** will be documented here.
  - Metabonaut now works with a stable version of Bioconductor (3.22)
- **Minor updates** can be found in the [News
  section](https://rformassspectrometry.github.io/Metabonaut/news/index.html).

------------------------------------------------------------------------

## 🎓 For R Beginners

The tutorials assume basic knowledge of **R** and **RMarkdown**. If
you’re new to these, we recommend starting with a short tutorial before
running the vignettes.

- **Learn Quarto (used for vignettes)**: [Quarto
  Guide](https://quarto.org/docs/guide/)
- **Learn RMarkdown**: [RMarkdown
  Book](https://bookdown.org/yihui/rmarkdown/)
- **Intro to R**: [Learn-R.org](https://learn-r.org/)
- **Interactive R course**:
  [Swirl](https://swirlstats.com/students.html)
- **Best Practices Cheatsheet**: [GitHub
  Repository](https://github.com/wurli/r-best-practice)

------------------------------------------------------------------------

## 🛠️ Known Issues

This is just the beginning of our **Metabonaut journey**, and we’re
actively refining the website. If you’re experiencing any issues:

✅ Ensure you have the latest versions of all required packages. 🐛 If
the issue persists, report it with a **reproducible example** on
[**GitHub
Issues**](https://github.com/rformassspectrometry/Metabonaut/issues).

Currently, there are **no known issues** with the code.

------------------------------------------------------------------------

## 🤝 Contribution

Interested in contributing? Please check out the [**RforMassSpectrometry
Contributions
Guide**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).

### 📜 Code of Conduct

We follow the [**RforMassSpectrometry Code of
Conduct**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct)
to maintain an inclusive and respectful community.

------------------------------------------------------------------------

## 🙌 Acknowledgements

![EU
Logo](https://github.com/rformassspectrometry/Metabonaut/raw/main/vignettes/images/EULogo.jpg)

EU Logo

This work is funded by the **European Union** under the
**HORIZON-MSCA-2021** project **101073062: HUMAN – Harmonising and
Unifying Blood Metabolic Analysis Networks**.

🔗 Learn more: [HUMAN Project Website](https://human-dn.eu/)
