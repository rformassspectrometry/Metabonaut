FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@outlook.com" \
      description="Docker container to run the different tutorials hosted on metabonaut. Includes Sirius 6.3 and RuSirius for the advanced feature annotation vignette. This version bases on the Bioconductor devel docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/
RUN rm -rf /home/rstudio/scripts /home/rstudio/vignettes/.quarto

## Global installation of required packages
## Need MsBackendMetaboLights to pre-download the dataset.
## Need MsIO because not on any repository
RUN Rscript -e "BiocManager::install(c('RforMassSpectrometry/MsIO', 'MsBackendMetaboLights', 'mzR') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

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

## Temporarily install SpectriPy from github
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/SpectriPy', ref = 'RELEASE_3_22', ask = FALSE, force = FALSE)"

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = c('Depends', 'Imports'), type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

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

