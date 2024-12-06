FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL name="rformassspectrometry/Metabonaut" \
      url="https://github.com/rformassspectrometry/Metabonaut" \
      maintainer="philippine.louail@eurac.edu" \
      description="Docker container to run the different tutorials hosted on metabonaut. This version bases on the Bioconductor devel docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Install quarto on docker image
RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

## Install the required packages
RUN Rscript -e "BiocManager::install('RforMassSpectrometry/MsIO', ask = FALSE, dependencies = TRUE)"

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.', dependencies = TRUE, type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"

RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} +
