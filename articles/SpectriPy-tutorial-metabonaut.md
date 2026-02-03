# LC-MS/MS Data Annotation using R and Python

## Introduction

The *SpectriPy* R package enables powerful mass spectrometry (MS) data
analysis workflows combining the strengths of Python and R MS libraries
([Graeve et al. 2025](#ref-graeve_spectripy_2025)). General concepts and
examples of the package are available in the package’s [*main*
vignette](https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html).

As an example for a combined R/Python workflow, we annotate the LC-MS/MS
spectra from the main end-to-end vignette using Python’s *matchms*
library.

The MS2 processing methods which will be demonstrated on the used
spectral reference library consists of the default filtering
functionality from *matchms* ([Huber et al.
2020](#ref-huber_matchms_2020)) (`default_filters()`,
`add_parent_mass()` and
[`normalize_intensities()`](https://rdrr.io/pkg/SpectriPy/man/filterSpectriPy.html)).

The MS2 spectral similarity algorithm which is demonstrated here is the
[`ModifiedCosine()`](https://rdrr.io/pkg/SpectriPy/man/compareSpectriPy.html)
from *matchms* ([Huber et al. 2020](#ref-huber_matchms_2020)).

The spectral reference library used for the annotation of the unknown
features in this tutorial originates from a small *in-house* reference
library (provided in MGF format) available in the *SpectriPy* package.

The analysis in this document is performed using R and Python code
chunks. A comment *\#’ R session:* or *\#’ Python session:* is used for
easier distinction of these.

## Load *SpectriPy*

Load the required R *SpectriPy* package. If you already have a Python
environment opened, please restart your Integrated Development
Environment and run this as a first command, to load the required
package *reticulate* and setup the conda environment ‘r-spectripy’
correctly. Please see the [Detailed information on installation and
configuration](https://rformassspectrometry.github.io/SpectriPy/articles/detailed-installation-configuration.html)
document for other options.

``` r

#' R session:

library(SpectriPy)
```

## Load query MS2 data

The LC-MS/MS query data used in this tutorial, are derived from the
Metabonaut resource ([Louail and Rainer
2025](#ref-louail_metabonaut_2025)). Introduction and a thorough
description of the preprocessing steps performed are described in [A
Complete End-to-End Workflow for untargeted LC-MS/MS Metabolomics Data
Analysis in
R](https://rformassspectrometry.github.io/Metabonaut/articles/end-to-end-untargeted-metabolomics.html).

First, we load the `Spectra` object with the MS2 spectra of the unknown
features found to be significant in the “Differential abundance
analysis”, see section [MS2-based
annotation](https://rformassspectrometry.github.io/Metabonaut/articles/end-to-end-untargeted-metabolomics.html#differential-abundance-analysis).
This `Spectra` object is shared as part of the *Metabonaut* package.

``` r

#' R session:

#' R MS package
library(Spectra)

#' Load the MS2 spectra of significant features
load(system.file("extdata", "spectra_significant_fts.RData",
                 package = "Metabonaut"))
ms2_ctr_fts
```

    MSn data (Spectra) with 315 spectra in a MsBackendMemory backend:
          msLevel     rtime scanIndex
        <integer> <numeric> <integer>
    1           2   147.357      2043
    2           2   148.587      2061
    3           2   149.817      2079
    4           2   152.297      2115
    5           2   147.376      2041
    ...       ...       ...       ...
    311         2   178.082      2481
    312         2   179.322      2499
    313         2   180.572      2517
    314         2   181.822      2535
    315         2   183.072      2553
     ... 39 more variables/columns.
    Processing:
     Filter: select retention time [10..240] on MS level(s) 1 2 [Tue Mar 18 11:56:42 2025]
     Filter: select MS level(s) 2 [Tue Mar 18 11:56:50 2025]
     Remove peaks based on their intensities and a user-provided function in spectra of MS level(s) 2. [Tue Mar 18 11:56:50 2025]
     ...19 more processings. Use 'processingLog' to list all. 

``` r

#' Print the available metadata, stored in the Spectra object
spectraVariables(ms2_ctr_fts)
```

     [1] "msLevel"                    "rtime"
     [3] "acquisitionNum"             "scanIndex"
     [5] "dataStorage"                "dataOrigin"
     [7] "centroided"                 "smoothed"
     [9] "polarity"                   "precScanNum"
    [11] "precursorMz"                "precursorIntensity"
    [13] "precursorCharge"            "collisionEnergy"
    [15] "isolationWindowLowerMz"     "isolationWindowTargetMz"
    [17] "isolationWindowUpperMz"     "peaksCount"
    [19] "totIonCurrent"              "basePeakMZ"
    [21] "basePeakIntensity"          "ionisationEnergy"
    [23] "lowMZ"                      "highMZ"
    [25] "mergedScan"                 "mergedResultScanNum"
    [27] "mergedResultStartScanNum"   "mergedResultEndScanNum"
    [29] "injectionTime"              "filterString"
    [31] "spectrumId"                 "ionMobilityDriftTime"
    [33] "scanWindowLowerLimit"       "scanWindowUpperLimit"
    [35] "electronBeamEnergy"         "mtbls_id"
    [37] "mtbls_assay_name"           "derived_spectral_data_file"
    [39] "collision_energy"           "feature_id"                

``` r

#' Print the feature_id of the first spectrum
ms2_ctr_fts$feature_id[1]
```

    [1] "FT0371"

## Filter query data

To ensure this `Spectra` object only contains MS2 data, we filter to
only MS2 spectra with more than 2 fragment peaks per spectrum using some
classical filtering functions from the *Spectra* package.

``` r

#' R session:

#' Filter MS2 level data
ms2_ctr_fts <- filterMsLevel(ms2_ctr_fts, 2L)

#' filter minimum 3 fragment peaks
ms2_ctr_fts <- ms2_ctr_fts[lengths(ms2_ctr_fts) >= 3]
ms2_ctr_fts
```

    MSn data (Spectra) with 291 spectra in a MsBackendMemory backend:
          msLevel     rtime scanIndex
        <integer> <numeric> <integer>
    1           2   147.357      2043
    2           2   148.587      2061
    3           2   152.297      2115
    4           2   147.376      2041
    5           2   148.616      2059
    ...       ...       ...       ...
    287         2   178.082      2481
    288         2   179.322      2499
    289         2   180.572      2517
    290         2   181.822      2535
    291         2   183.072      2553
     ... 39 more variables/columns.
    Processing:
     Filter: select retention time [10..240] on MS level(s) 1 2 [Tue Mar 18 11:56:42 2025]
     Filter: select MS level(s) 2 [Tue Mar 18 11:56:50 2025]
     Remove peaks based on their intensities and a user-provided function in spectra of MS level(s) 2. [Tue Mar 18 11:56:50 2025]
     ...20 more processings. Use 'processingLog' to list all. 

## Load reference MS2 data

As the *in-house* spectral library, we import a small test data file in
MGF format. This file is part of the *SpectriPy* package and we below
define its path on the local computer.

``` r

#' R session:

#' Define a variable with the path and file name of the MGF file
mgf_file <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
```

The loading of this file is then performed using the Python *matchms*
library. The variable with the MGF file name can be accessed in the
associated Python session using `r.mgf_file`. The loaded object is a
Python list of `matchms.Spectrum` objects.

``` python
#' Python session:

from matchms.importing import load_from_mgf

#' Read spectra from an MGF formatted file, as Spectrum object
mgf_py = list(load_from_mgf(r.mgf_file))

#' Nr of spectra
len(mgf_py)
```

    100

``` python
#' Access the first spectrum
mgf_py[0]
```

    Spectrum(precursor m/z=259.06, 3 fragments between 213.1 and 259.1)

Note that we can also access the first spectrum from an R session, by
starting the command with `py$`.

``` r

#' R session:

#' Access the first spectrum
py$mgf_py[[1]]
```

    Spectrum(precursor m/z=259.06, 3 fragments between 213.1 and 259.1)

## Translate query MS data to Python

First, we check if the R `Spectra` object containing the query MS2 data
can be accessed in Python using the `r.` prefix.

``` python
#' Python session:

#' check if the r Spectra object can be accessed in python using the 'r.'
#' prefix. Print the first spectrum
r.ms2_ctr_fts[1]
```

    Spectrum(precursor m/z=138.05, 4 fragments between 73.1 and 92.1)

``` python
#' Show which metadata is available in the first spectrum
r.ms2_ctr_fts[0].metadata.keys()
```

    dict_keys(['precursor_mz', 'precursor_intensity', 'charge', 'retention_time', 'collision_energy', 'isolation_window_target_mz', 'ms_level'])

Second, we translate the `Spectra` object `ms2_ctr_fts` to respective
Python `matchms.Spectrum` objects using the
[`rspec_to_pyspec()`](https://rdrr.io/pkg/SpectriPy/man/conversion.html)
function. With the
[`py_set_attr()`](https://rstudio.github.io/reticulate/reference/py_set_attr.html)
function we assign the variable directly to an attribute in the Python
session (which avoids repeated cross-programming language references).

``` r

#' R session:

#' Add mapping for additional spectra variables to the default mapping in R and
#' python, respectively
map = c(defaultSpectraVariableMapping(),
        feature_id = 'feature_id')

#' Convert to py Spectrum
py_set_attr(py, "ms2_ctr_fts_py", rspec_to_pyspec(ms2_ctr_fts, mapping = map))
```

We can now access the converted `Spectra` object in Python.

``` python
#' Python session:

#' Access the first converted spectrum
ms2_ctr_fts_py[0]
```

    Spectrum(precursor m/z=138.05, 4 fragments between 73.1 and 92.1)

## Filter the reference library

Before we run the spectral comparisons of our query data to the MGF
reference library, we first apply some MS2 processing from *matchms*.
Default filtering from *matchms* is performed to standardize ion mode,
correct charge and more. See the [matchms filtering
documentation](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html).
Also, we keep only reference spectra with a precursor m/z as the
similarity algorithm we use later requires spectra to have a precursor
m/z.

``` python
#' Python session:

from matchms.filtering import default_filters, normalize_intensities, add_parent_mass

#' Apply filters to clean and enhance each spectrum
clean_mgf_py = []
for spectrum in mgf_py:
    #' Apply default filter to standardize ion mode, correct charge and more.
    #' Default filter is fully explained at
    #' https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html
    spectrum = default_filters(spectrum)
    #' For missing precursor_mz field: check if there is “pepmass”” entry instead
    spectrum = add_parent_mass(spectrum)
    #' Scale peak intensities to maximum of 1
    spectrum = normalize_intensities(spectrum)
    #' Only add spectra that have a precursor m/z
    if "precursor_mz" in spectrum.metadata:
        clean_mgf_py.append(spectrum)

#' Nr of spectra
len(clean_mgf_py)
```

    78

## Calculating spectra similarities using the Modified Cosine algorithm from *matchms*

We calculate the pairwise spectral similarity between the query spectra
and the reference library spectra using Python’s *matchms* library.

Here, we use the spectral similarity algorithm ModifiedCosine from
*matchms*, from source
[matchms](https://github.com/matchms/matchms/blob/master/README.rst).
This algorithm can easily be exchanged for another spectral similarity
calculation from *matchms* (see
[here](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html)).

``` python
#' Python session:

from matchms import calculate_scores
from matchms.similarity import ModifiedCosine

#' Calculate Cosine similarity scores between all spectra
#' For other similarity score methods see
#' https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html
similarity_score = ModifiedCosine(tolerance = 0.1)
scores = calculate_scores(references = clean_mgf_py,
                          queries = ms2_ctr_fts_py,
                          similarity_function = similarity_score)
```


    Calculating similarities:   0%|          | 0/78 [00:00<?, ?it/s]
    Calculating similarities:   1%|1         | 1/78 [00:01<02:28,  1.93s/it]
    Calculating similarities:  12%|#1        | 9/78 [00:02<00:11,  5.94it/s]
    Calculating similarities:  21%|##        | 16/78 [00:02<00:05, 11.63it/s]
    Calculating similarities:  31%|###       | 24/78 [00:02<00:02, 19.22it/s]
    Calculating similarities:  41%|####1     | 32/78 [00:02<00:01, 27.57it/s]
    Calculating similarities:  51%|#####1    | 40/78 [00:02<00:01, 35.62it/s]
    Calculating similarities:  62%|######1   | 48/78 [00:02<00:00, 41.76it/s]
    Calculating similarities:  71%|#######   | 55/78 [00:02<00:00, 47.40it/s]
    Calculating similarities:  81%|########  | 63/78 [00:02<00:00, 53.68it/s]
    Calculating similarities:  91%|#########1| 71/78 [00:02<00:00, 59.69it/s]
    Calculating similarities: 100%|##########| 78/78 [00:02<00:00, 26.06it/s]

``` python
scores
```

    <78x291x2 stacked sparse array containing scores for ('ModifiedCosine_score', 'ModifiedCosine_matches') with 12376 stored elements in COOrdinate format>

We next rearrange the spectra similarity results to make a data frame
containing the best matched compound name (derived from the reference
library) per queried spectrum.

First, we extract and transpose the scores as a python array. Each row
of the array will contain the similarity scores of one spectrum from our
query spectra `ms2_ctr_fts_py` against the cleaned reference library
`clean_mgf_py`.

``` python
#' Python session:

#' Convert to array and transpose
sim_matchms = scores.to_array()["ModifiedCosine_score"]
sim_matchms = sim_matchms.T

#' Contains 1 row for each spectrum in query
sim_matchms.shape
```

    (291, 78)

Next, we create a data frame with per queried spectrum from our unknown
variables, the compound name of the higest matching spectra from the
reference library and the corresponding similarity score.

``` python
#' Python session:

import numpy as np
import pandas as pd

#' Prepare results list
results = []
for i in range(sim_matchms.shape[0]):
    #' row is the query, keep nr in the results instead of replacing by eg id
    name_row = ms2_ctr_fts_py[i].get('feature_id')
    row_values = sim_matchms[i].copy()
    #' match with higest col nr from the references
    max_col = np.argmax(row_values)  #' Find column index of max value
    max_value = row_values[max_col]  #' Get max value
    #' replace the nr of refererences with the name
    name_max_col = clean_mgf_py[max_col].get('compound_name')
    results.append({"query": i + 1,
                    "query_feature_id": name_row,
                    "reference": max_col,
                    "reference_compound_name": name_max_col,
                    "ModifiedCosine_score": max_value})


#' Convert to DataFrame
df = pd.DataFrame(results)

#' Print the first 5 rows of the unfiltered DataFrame
df.head()
```

       query query_feature_id  ...  reference_compound_name ModifiedCosine_score
    0      1           FT0371  ...         Benzyl succinate             0.554071
    1      2           FT0371  ...         Benzyl succinate             0.557885
    2      3           FT0371  ...      L-(+)-Ergothioneine             0.103269
    3      4           FT0371  ...         Benzyl succinate             0.464949
    4      5           FT0371  ...         Benzyl succinate             0.508411

    [5 rows x 5 columns]

\[!\] **Caution**: As the higest score is taken as criteria for the
annotation, a lot of caution is needed evaluation the trueness of the
match. A low score is not reliable, as the similarity algorithm will
calculate a score for each pair of spectra. Therefore, a match will
always be found. In addition, if your unknown compound is absent in the
reference library, it will match wrongly to another compound that is
present in the database.

Below we filter the results retaining only matches with a similarity
value above 0.6. Above this value, the potential annotations need to be
validated using e.g. rerunning samples in the presence of commercial
standards.

``` python
#' Python session:

#' Keep only rows where score >= 0.6
df_filtered = df[df["ModifiedCosine_score"] >= 0.6]
```

``` r

#' R session:

library(pander)
#' Print the filtered DataFrame
pandoc.table(py$df_filtered, style = "rmarkdown", split.table = Inf)
```

|   | query | query_feature_id | reference | reference_compound_name | ModifiedCosine_score |
|:--:|:--:|:--:|:--:|:--:|:--:|
| **14** | 15 | FT0371 | 52 | Nordihydroguaiaretic acid | 0.6667 |
| **33** | 34 | FT0565 | 35 | Ethambutol | 0.9623 |
| **34** | 35 | FT0565 | 35 | Ethambutol | 0.8809 |
| **37** | 38 | FT0565 | 35 | Ethambutol | 0.9412 |
| **38** | 39 | FT0565 | 35 | Ethambutol | 0.8681 |
| **40** | 41 | FT0565 | 35 | Ethambutol | 0.935 |
| **41** | 42 | FT0565 | 35 | Ethambutol | 0.7319 |
| **44** | 45 | FT0565 | 35 | Ethambutol | 0.89 |
| **48** | 49 | FT0565 | 35 | Ethambutol | 0.6822 |
| **51** | 52 | FT0565 | 38 | Benzyl succinate | 0.7218 |
| **52** | 53 | FT0565 | 38 | Benzyl succinate | 0.7516 |
| **53** | 54 | FT0565 | 47 | Bisphenol AF | 0.6977 |
| **56** | 57 | FT0565 | 39 | Isethionate | 0.6666 |
| **57** | 58 | FT0565 | 38 | Benzyl succinate | 0.7443 |
| **67** | 68 | FT0732 | 72 | Prednisolone_Tebutate | 0.7048 |
| **68** | 69 | FT0732 | 72 | Prednisolone_Tebutate | 0.7729 |
| **70** | 71 | FT0732 | 72 | Prednisolone_Tebutate | 0.7103 |
| **71** | 72 | FT0732 | 72 | Prednisolone_Tebutate | 0.6087 |
| **73** | 74 | FT0732 | 72 | Prednisolone_Tebutate | 0.6523 |
| **74** | 75 | FT0732 | 72 | Prednisolone_Tebutate | 0.8334 |
| **87** | 88 | FT0732 | 19 | Isoproturon-didemethyl | 0.6235 |
| **89** | 90 | FT0732 | 72 | Prednisolone_Tebutate | 0.6793 |
| **97** | 98 | FT0732 | 9 | N-Methyl-2-pyrrolidone | 0.7953 |
| **103** | 104 | FT0732 | 72 | Prednisolone_Tebutate | 0.6041 |
| **104** | 105 | FT0732 | 72 | Prednisolone_Tebutate | 0.6836 |
| **109** | 110 | FT0732 | 72 | Prednisolone_Tebutate | 0.6295 |
| **111** | 112 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6005 |
| **112** | 113 | FT0732 | 16 | Caffeine | 0.6181 |
| **113** | 114 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6843 |
| **125** | 126 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6863 |
| **126** | 127 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7172 |
| **127** | 128 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7157 |
| **128** | 129 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.721 |
| **129** | 130 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6315 |
| **134** | 135 | FT0732 | 72 | Prednisolone_Tebutate | 0.7612 |
| **137** | 138 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6478 |
| **139** | 140 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6301 |
| **140** | 141 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7013 |
| **141** | 142 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7017 |
| **142** | 143 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7018 |
| **143** | 144 | FT0732 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7863 |
| **186** | 187 | FT0732 | 9 | N-Methyl-2-pyrrolidone | 0.636 |
| **220** | 221 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.737 |
| **221** | 222 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.733 |
| **222** | 223 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7256 |
| **223** | 224 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7305 |
| **227** | 228 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.6283 |
| **228** | 229 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7359 |
| **229** | 230 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7296 |
| **230** | 231 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.73 |
| **239** | 240 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7352 |
| **240** | 241 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7351 |
| **247** | 248 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7391 |
| **248** | 249 | FT0845 | 11 | 3-(4-Isopropylphenyl)isobutyraldehyde | 0.7332 |

To visually inspect how good the query and reference spectra match, we
refer to the [Metabonaut
resource](https://rformassspectrometry.github.io/Metabonaut/articles/end-to-end-untargeted-metabolomics.html)
on how to generate the mirror plots and perform precursor *m/z*
filtering (e.g. maximum 1 Da difference).

## Session information

``` r

#' R session:

sessionInfo()
```

    R version 4.5.2 (2025-10-31)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.3 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

    time zone: Etc/UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods
    [8] base

    other attached packages:
    [1] pander_0.6.6        Spectra_1.20.1      BiocParallel_1.44.0
    [4] S4Vectors_0.48.0    BiocGenerics_0.56.0 generics_0.1.4
    [7] SpectriPy_1.0.0     reticulate_1.44.1

    loaded via a namespace (and not attached):
     [1] cli_3.6.5              knitr_1.51             rlang_1.1.7
     [4] xfun_0.56              ProtGenerics_1.42.0    otel_0.2.0
     [7] png_0.1-8              jsonlite_2.0.0         clue_0.3-66
    [10] rprojroot_2.1.1        htmltools_0.5.9        rmarkdown_2.30
    [13] grid_4.5.2             evaluate_1.0.5         MASS_7.3-65
    [16] fastmap_1.2.0          yaml_2.3.12            IRanges_2.44.0
    [19] MsCoreUtils_1.22.1     cluster_2.1.8.1        compiler_4.5.2
    [22] codetools_0.2-20       fs_1.6.6               Rcpp_1.1.1
    [25] here_1.0.2             MetaboCoreUtils_1.18.1 lattice_0.22-7
    [28] digest_0.6.39          parallel_4.5.2         Matrix_1.7-4
    [31] withr_3.0.2            tools_4.5.2           

## References

Graeve, Marilyn De, Wout Bittremieux, Thomas Naake, Carolin Huber,
Matthias Anagho-Mattanovich, Nils Hoffmann, Pierre Marchal, et al. 2025.
“SpectriPy: Enhancing Cross-Language Mass Spectrometry Data Analysis
with R and Python.” *Journal of Open Source Software* 10 (109): 8070.
<https://doi.org/10.21105/joss.08070>.

Huber, Florian, Stefan Verhoeven, Christiaan Meijer, Hanno Spreeuw,
Efraín Castilla, Cunliang Geng, Justin Van Der Hooft, et al. 2020.
“Matchms - Processing and Similarity Evaluation of Mass Spectrometry
Data.” *Journal of Open Source Software* 5 (52): 2411.
<https://doi.org/10.21105/joss.02411>.

Louail, Philippine, and Johannes Rainer. 2025.
“Rformassspectrometry/Metabonaut: Metabonaut Version 1.0.0.” Zenodo.
<https://doi.org/10.5281/ZENODO.15062930>.
