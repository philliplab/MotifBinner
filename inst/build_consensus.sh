#!/bin/bash

R -e "library(knitr); knit2html('consensus_investigation.Rmd')"

#firefox consensus_investigation.html
