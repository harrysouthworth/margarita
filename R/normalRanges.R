#' Normal ranges for some common laboratory variables
#' 
#' @name normalRanges
#' @docType data
#' @details The normal ranges for some common laboratory variables, as obtained
#'   from the Merck Manual, January 2015. In each case, 'lln' refers to the lower
#'   limit of normal and 'uln' to the upper limit of normal.
#' @author Harry Southworth
#' @references Merck Manual, http://www.merckmanuals.com/professional/appendixes/normal_laboratory_values/blood_tests_normal_values.html#v8508814

if (FALSE){
  test <- c("alanine aminotransferase", "aspartate aminotransferase", "bilirubin, total",
            "alkaline phosphatase", "calcium", "creatinine", "fasting glucose", "albumin",
            "hemoglobin a1c", "potassium", "protein, total", "sodium", "urea nitrogen", "uric acid")
  
  abbreviation <- c("ALT", "AST", "TBL", "ALP", "Calcium", "Creat", "Gluc", "Alb", "HbA1c", "K",
                    "Protein", "Na", "BUN", "UA")
  specimen <- c(rep("serum", 6), "plasma", "serum", "blood", "serum", "serum", "serum", "serum", "serum")
  
  lln.conventional <- c(0, 0, 0.3, 36, 9, .7, 70, 3.5, 4.7, 3.5, 6, 136, 8, 2.5)
  uln.conventional <- c(35, 35, 1.2, 92, 10.5, 1.3, 105, 5.5, 8.5, 5, 7.8, 145, 20, 8)
  
  lln.si <- c(0, 0, 5.1, .5, 2.2, 61.9, 3.9, 35, 4.7, 3.5, 60, 136, 2.9, .15)
  uln.si <- c(.58, .58, 20.5, 1.5, 2.6, 115, 5.8, 55, 8.5, 5, 78, 145, 7.1, .47)
  
  conventional.units <- c("U/L", "U/L", "mg/dL", "U/L", "mg/dL", "mg/dL", "mg/dL", "g/dL", "%", "mEq/L", "g/dL", "mEq/L", "mg/dL", "mmol/L")
  si.units <- c("ukat/L", "ukat/L", "umol/L", "ukat/L", "mmol/L", "umol/L", "mmol/L", "g/L", "%", "mmol/L", "g/L", "mmol/L", "mg/dL", "mmol/L")
  
  normalRanges <- data.frame(test=test, abbreviation=abbreviation, specimen=specimen,
                             lln.conventional=lln.conventional, uln.conventional=uln.conventional,
                             lln.si=lln.si, uln.si=uln.si,
                             conventional.units=conventional.units, si.units=si.units)
  save(normalRanges, file="~/Work/repos/github/margarita/data/normalRanges.rda")
}
