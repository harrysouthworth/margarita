suppressMessages(library(margarita))
suppressMessages(library(xtable))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

options(width=60)

seed <- as.numeric(gsub("-", "", Sys.Date())) # no need to use date here
set.seed(seed)
