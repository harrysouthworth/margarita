# Source:
# http://www.merckmanuals.com/professional/appendixes/normal_laboratory_values/blood_tests_normal_values.html#v8508814
# at 2014-09-30

test <- c("alanine aminotransferase", "aspartate aminotransferase", "bilirubin, total",
          "alkaline phosphatase")
abbreviation <- c("ALT", "AST", "TBL", "ALP")
specimen <- rep("serum", length(test))

lln.conventional <- c(0, 0, 0.3, 36)
uln.conventional <- c(35, 35, 1.2, 92)

lln.si <- c(0, 0, 5.1, .5)
uln.si <- c(.58, .58, 20.5, 1.5)

conventional.units <- c("U/L", "U/L", "mg/dL", "U/L")
si.units <- c("ukat/L", "ukat/L", "umol/L", "ukat/L")

normalRanges <- data.frame(test=test, abbreviation=abbreviation, specimen=specimen,
                           lln.conventional=lln.conventional, uln.conventional=uln.conventional,
                           lln.si=lln.si, uln.si=uln.si,
                           conventional.units=conventional.units, si.units=si.units)