#! /bin/bash

cd /home/harry/Work/repos/github/margarita/vignettes
mv margarita.knitr margarita.Rnw
mv safety.knitr safety.Rnw

R < build.R --no-save

mv margarita.Rnw margarita.knitr
mv safety.Rnw safety.knitr

make all
