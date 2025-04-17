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

## Use SpectriPy with virtual env to avoid need to install miniconda
ENV SPECTRIPY_USE_CONDA="FALSE"

## Install SpectriPy and caching files for rstudio user
USER rstudio
RUN Rscript -e "install.packages('reticulate')" && \
    Rscript -e "BiocManager::install('RforMassSpectrometry/SpectriPy')" && \
    Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())"

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = TRUE, type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())" && \
    find vignettes/ -name "*.html" -type f -delete && \
    find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*

## root user needed for rstudio server properly working
USER root
