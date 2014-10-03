if (dataRoot == ""){ # Run with texmex liver data
  lab <- lb
  dem <- dm  
} else {
  lab <- read.sas7bdat(paste0(dataPath, "/", labs))
  dem <- read.sas7bdat(paste0(dataPath, "/", demog))
}

# Get number of unique subjects
lab <- lab[lab[, popflag] %in% popyes, ]
dem <- dem[dem[, popflag] %in% popyes, ]

nSub <- length(!is.na(unique(lab[, subject])))

# Subset to lab variable of interest
lab <- lab[lab[, test] == testval, ]

# Add arm if it's missing
lab <- addVariables(lab, dem, subject=subject, vars=arm)

# Add baseline if it's missing
lab <- getBaselines(lab, subject=subject, visit=visit, baseline.visit=basevisit,
                    value=value, baseline=baseline)

# Get maxima, or minima. XXX Need to tidy this
if (minmax == "max"){
  aggfun <- max
} else if (minmax == "min"){
  aggrun <- min
} else{
  stop("aggregation function must be either min or max")
}
lab <- getAggregateData(lab, subject=subject, visit=visit,
                        baseline.visit=basevisit, max.study.visit=maxvisit,
                        baseline=baseline, value=value, aggregate.fun=aggfun)

if (nrow(lab) != length(unique(lab[, subject])))
  stop("Number of rows in aggregated data does not equal number of study subjects")

if (minmax == "min") lab[, value] <- -lab[, value]

# Add treatment group information
lab <- addVariables(lab, dem, subject=subject, vars=arm)

# Catch and remove missing values
# XXX Need a function to do this neatly
lab <- lab[!(is.na(lab[, arm]) | is.na(lab[, subject]) | is.na(lab[, baseline]) | is.na(lab[, value])), ]

nMiss <- nSub - length(unique(lab[, subject]))
if (nMiss == 0) nMiss <- "no"
