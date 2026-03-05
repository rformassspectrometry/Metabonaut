#!/usr/bin/with-contenv bash
# shellcheck shell=bash
WD="/tmp/.sirius"

mkdir $WD
/home/rstudio/sirius/bin/sirius --workspace=$WD service --port=9999 &

chmod -R a+rw $WD
