FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@outlook.com" \
      description="Docker container to run the different tutorials hosted on metabonaut. Includes Sirius 6.3 and RuSirius for the advanced feature annotation vignette. This version bases on the Bioconductor devel docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

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

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = c('Depends', 'Imports'), type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

## root user needed for rstudio server properly working
USER root

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*

## Install Miniconda (needed for Sirius 6.3 installation for interactive use)
RUN apt-get update && apt-get install -y curl && \
    curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda clean --all -y && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile && \
    echo "conda activate base" >> ~/.bashrc

## Configure conda channels
RUN /opt/conda/bin/conda config --remove channels defaults || true && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --set channel_priority strict

## Install Sirius 6.3 via conda (required by RuSirius)
RUN /opt/conda/bin/conda install --override-channels -c conda-forge r-sirius-ms -y

## Set Sirius on the PATH so R/RuSirius can find the executable
RUN SIRIUS_PATH=$(find /opt/conda/pkgs -maxdepth 1 -type d -name "sirius-ms-*" | head -1) && \
    echo "export PATH=${SIRIUS_PATH}/bin:\${PATH}" >> /etc/profile.d/sirius.sh && \
    echo "export PATH=${SIRIUS_PATH}/bin:\${PATH}" >> /home/rstudio/.bashrc

## Install RuSirius (R interface to Sirius) for interactive use
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/RuSirius', ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"
