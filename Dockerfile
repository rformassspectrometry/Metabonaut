FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@outlook.com" \
      description="Docker container to run the different tutorials hosted on metabonaut. Includes Sirius 6.3 and RuSirius for the advanced feature annotation vignette. This version is based on the Bioconductor release 3.23 docker image." \
      license="CC-BY-SA-4.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/
RUN rm -rf /home/rstudio/scripts /home/rstudio/vignettes/.quarto

## Global installation of required packages
## Need MsBackendMetaboLights to pre-download the dataset.
## Need MsIO because not on any repository
RUN Rscript -e "install.packages('remotes'); BiocManager::install(c('RforMassSpectrometry/MsIO', 'MsBackendMetaboLights', 'mzR') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

## Use SpectriPy with virtual env to avoid need to install miniconda
ENV SPECTRIPY_USE_CONDA="FALSE"

## Install SpectriPy and caching files for rstudio user
## NOTE: installing SpectriPy from GitHub (devel) for now instead of the stable
## Bioconductor release; revert to 'SpectriPy' once the Bioc version is fixed.
USER rstudio
RUN Rscript -e "install.packages('reticulate')" && \
    Rscript -e "BiocManager::install('RforMassSpectrometry/SpectriPy', ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)" && \
    Rscript -e "library(MsBackendMetaboLights);Spectra('MTBLS8735', source = MsBackendMetaboLights())"

## Install the current package and build its vignettes in two steps.
## Step 1: install Metabonaut + all Suggests WITHOUT building any vignettes.
##   dependencies = TRUE is needed because the vignettes use Suggests packages;
##   build_vignettes = FALSE avoids building the *dependencies'* own vignettes
##   (e.g. RuSirius's vignettes require a running Sirius instance and would fail).
RUN Rscript -e "remotes::install_local('.', dependencies = TRUE, type = 'source', build_vignettes = FALSE, repos = BiocManager::repositories())"
## Step 2: rebuild ONLY Metabonaut with its vignettes (deps already installed).
##   force = TRUE is required, otherwise remotes skips the already-installed same
##   version and never builds the vignettes. chmod makes the read-only alabaster
##   objects in inst/extdata writable so R CMD build can clean up its temp copy.
RUN chmod -R u+w . && \
    Rscript -e "remotes::install_local('.', dependencies = FALSE, type = 'source', build_vignettes = TRUE, force = TRUE, repos = BiocManager::repositories())"

## Fail the build if the vignettes were not actually built and installed.
RUN Rscript -e "vi <- tools::getVignetteInfo('Metabonaut'); if (nrow(vi) < 1L) stop('No vignettes were installed for Metabonaut - vignette build failed'); message('OK: ', nrow(vi), ' vignette(s) installed:'); print(unname(vi[, 'File']))"

## root user needed for rstudio server properly working
USER root

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*

## Install sirius.
RUN wget -nv https://github.com/sirius-ms/sirius/releases/download/v6.3.3/sirius-6.3.3-linux-x64.zip && \
    unzip sirius-*.zip && \
    rm sirius-*.zip && \
    chown -R rstudio:rstudio sirius && \
    ln -s /home/rstudio/sirius/bin/sirius /usr/local/bin/sirius && \
    echo "export PATH=/home/rstudio/sirius/bin:$PATH" >> /home/rstudio/.bashrc

COPY ./scripts/sirius-init.sh /etc/cont-init.d/03_sirius
RUN chmod a+x /etc/cont-init.d/03_sirius

## Install RuSirius (R interface to Sirius) for interactive use
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/RuSirius', ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

