FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@outlook.com" \
      description="Docker container to run the different tutorials hosted on metabonaut. This version bases on the Bioconductor devel docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Global installation of required packages
## Need MsBackendMetaboLights to pre-download the dataset.
## Need MsIO because not on any repository
RUN Rscript -e "BiocManager::install(c('RforMassSpectrometry/MsIO', 'MsBackendMetaboLights', 'mzR', 'RforMassSpectrometry/MsBackendMassbank') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

## Install keyring package from github using pak
RUN Rscript -e "install.packages('pak');pak::pak('r-lib/keyring', ask = FALSE)"

## Use SpectriPy with virtual env to avoid need to install miniconda
ENV SPECTRIPY_USE_CONDA="FALSE"

## Install SpectriPy and caching files for rstudio user
USER rstudio
#RUN Rscript -e "install.packages('reticulate')" && \
#    Rscript -e "BiocManager::install('RforMassSpectrometry/SpectriPy', dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)" && \
#    Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())"

RUN Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())"

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = c('Depends', 'Imports'), type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

## root user needed for rstudio server properly working
USER root

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*
