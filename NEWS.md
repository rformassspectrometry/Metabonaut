# Metabonaut 1.2

## Changes in 1.2.0

- Addition of a new vignette depicting a large scale processing using xcms.
  Due to its size, this is added as a pre-computed vignette to prevent re-knitting
  as much as possible.


# Metabonaut 1.1

## Changes in 1.1.2

- Adding authors that have contributed to specific method/vignettes to the
  respective vignette's header.
- Update README.md to add the link to SpectriPy vignette.
- fix ORCID  in description so it gets recognize as an ORCID.

## Changes in 1.1.1

- Update docker to skip installation of R packages and caching of data files for
  the *root* user.
- Fix GHA to properly install all Python libraries during `library(SpectriPy)`.

## Changes in 1.1.0

- Add a vignette for combined R/Python data analysis using *reticulate* and
  *SpectriPy*.

## Changes in 1.1.0

- Addition of a devel branch to the repository to maintain a stable version.

## Changes in 1.0.1

- Minor fixes of DESCRIPTION and .yml
- Addition of a devel branch to the repository to maintain a stable version.

## Changes in 1.0.0

- More plots in the data investigation vignette.
- In the main end to end vignette: added saving of the Spectra object with the
  significant features associated MS1 and MS2 data. To be used later by other
  vignettes.


# Metabonaut 0.0

## Changes in 0.0.7

- Fix typo and cut-off for MS2 annotation.
- Updated figure numbering.
- Added acknowledgment to community.
- added a figure in retention time alignment vignette.

## Changes in 0.0.6

- Moving The PercentMissing filtering at the end of
  the pre-processing steps as it needs to be done
  before normalization.
- Addition  of collapsing code to improve readability.
- Reduction of table size

## Changes in 0.0.5

- Require *MsIO* version 0.0.8 to allow reading of stored
  `MsBackendMetaboLights` objects.
- Small updates and changes to the Seamless Alignment vignette:
  - simplify import of previously stored result object
  - avoid using the variable name `param` for every parameter object

## Changes in 0.0.4
- Required *alabaster.se*.
- In the end-to-end vignette:
  - Removal of Spectra data in depth visualisation to
    move to the Data investigation vignette
  - Removal of internal standard matching to features
    in the Normalization part.
- Save an load *lcms1* and *res* object from the end-to-end workflow
  to be used in the Seamless Alignment vignette. Using *MsIO* and
  *alabaster.se*

## Changes in 0.0.3
- Require *xcms* version 4.3.4 and install the package from github.

## Changes in 0.0.2
- Switch to Quarto instead of Rmarkdown
- Addition of Alignment to reference dataset vignette
- Addition of the Data investigation vignette
- Addition of the Install vignette

## Changes in 0.0.1
- Addition of basic files for a workflow package.
- Addition of the end-to-end vignette.
