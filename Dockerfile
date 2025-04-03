FROM bioconductor/bioconductor_docker:devel
## update to 3.21 in april to have stable versioning.

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@eurac.edu" \
      description="Docker container to run the different tutorials hosted on metabonaut. This version bases on the Bioconductor devel docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Global installation of required packages
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/MsBackendMetaboLights', ask = FALSE, dependencies = TRUE)" && \
    Rscript -e "BiocManager::install('RforMassSpectrometry/MsIO', ask = FALSE, dependencies = TRUE)"

## This below would actually work.
## Create the BiocFileCache and cache the data files to avoid repeated downloads
# USER rstudio
# RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
#     && bash Miniforge3.sh -b -p "/home/rstudio/conda" \
#     && rm Miniforge3.sh \
#     && . "/home/rstudio/conda/etc/profile.d/conda.sh" \
#     && conda init \
#     && Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())" \
#     && Rscript -e "install.packages('reticulate')"
# USER root

ENV SPECTRIPY_USE_CONDA="FALSE"
## alternative echo "Sys.setenv(SPECTRIPY_USE_CONDA='FALSE')" > ${HOME}/.Rprofile && \

## Alternative installation using Python envs
USER rstudio
RUN Rscript -e "install.packages('reticulate')" && \
    Rscript -e "BiocManager::install('RforMassSpectrometry/SpectriPy')" && \
    Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())"

## USER root

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = TRUE, type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} +
