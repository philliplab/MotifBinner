#!/bin/bash

R -e "library(knitr); knit2html('mislabel_investigation.Rmd')"

firefox mislabel_investigation.html
