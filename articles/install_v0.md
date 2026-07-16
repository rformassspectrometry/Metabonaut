# Install

For manual installation, an [R](https://www.r-project.org/) version \>=
4.4.0 is required.

## Running workflows in the cloud

No installation is needed to run the workflows: Metabonaut is available
as a *Curated Workshop* on the [Bioconductor Workshop
platform](https://workshop.bioconductor.org), a Galaxy-based service
hosted on the Jetstream2 academic cloud (NSF-ACCESS).

- Log in on
  [workshop.bioconductor.org](https://workshop.bioconductor.org) with a
  Bioconductor account.
- Launch the **Metabonaut** workshop. This starts a private RStudio
  session, based on the same docker image as below, with all required
  packages and the example data already installed.
- Open any of the Quarto files in the *vignettes* folder and evaluate
  the R code blocks.

Note that these sessions are ephemeral: all data and progress are erased
when the workshop is stopped, so download a local copy of anything you
want to keep.

## Running workflows locally

To install on your computer all the packages necessary for the workflows
run the code as follow:

``` r

install.packages("BiocManager")
BiocManager::install('RforMassSpectrometry/MsIO', ask = FALSE,
                     dependencies = TRUE)

BiocManager::install("RforMassSpectrometry/Metabonaut",
                     dependencies = TRUE, ask = FALSE, update = TRUE)
```

If you get some error message when downloading form GitHub it might be
due to an expired token. Remove them using this code below and try to
install again:

``` r

#See if any token is there:
gitcreds::gitcreds_get()

# Delete them:
gitcreds::gitcreds_delete()
```

## Docker image

The vignettes files along with an R runtime environment including all
required packages and the RStudio (Posit) editor are all bundled in a
*docker* container.

After installation, this docker container can be run on the computer and
the code and examples from the vignettes can be evaluated within this
environment (without the need to install any additional packages or
files).

- If you don’t already have, install [docker](https://www.docker.com/).
  Find installation information
  [here](https://docs.docker.com/desktop/).
- Get the [docker
  image](https://hub.docker.com/r/rformassspectrometry/metabonaut) of
  this tutorial e.g. from the command line with:

&nbsp;

    docker pull rformassspectrometry/metabonaut:latest

- Start the docker container, either through the Docker Desktop, or on
  the command line with

&nbsp;

    docker run -e PASSWORD=bioc -p 8787:8787 rformassspectrometry/metabonaut:latest

- Enter [`http://localhost:8787`](http://localhost:8787) in a web
  browser and log in with username `rstudio` and password `bioc`.
- In the RStudio server version: open any of the Quarto files in the
  *vignettes* folder and evaluate the R code blocks in that document.
