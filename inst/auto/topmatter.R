
######################### Information on data location #########################
drug <- "D1234" # used to construct path to data
study <- "0001" # used to construct path to data
dataRoot <- "" # default to the root of the directory tree containing the data
labs <- "lb.sas7bdat"
demog <- "dm.sas7bdat"

####################### Information on SAS data contents #######################
popflag <- "saffl"
popyes <- c("Y", "YES", "Yes", "y", "yes")
subject <- "usubjid"
visit <- "visitnum"
baseline <- "base"
basevisit <- 1
maxvisit <- Inf # values of visit > maxvisit will be dropped
test <- "paramcd"
testval <- "ALT"
value <- "aval"
minmax <- "max" # character string, anticipating a user interface
arm <- "arm" # could be a vector of length > 1, e.g. drug, dose

############################## Information for EVM #############################
trans <- "log" # character string, anticipating a user interface
models <- list(paste("~", arm)) # character string, anticipating a user interface
returnLevels <- c(100, 500, 1000)
ULN <- 35 # can be dropdown from normalRanges object or something
multiplesOfULN <- c(3, 5, 10, 20) # defaults can be inferred from CTC object or something

# XXX Next 3 are for speed whilst in dev. Get rid and use good defaults.
thin <- 1
burn <- 500
iter <- 10500 

