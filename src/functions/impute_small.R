impute_small <- function(dataframe,number_of_values) {
  tmp_impute <- min(dataframe)
  return(rnorm(number_of_values, tmp_impute/2, tmp_impute/8))
}