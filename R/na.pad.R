na.pad <- function (vec, ind) {
  vec0 <- numeric(length(ind))
  vec0[!ind] <- vec
  vec0[ind] <- NA
  return(vec0)
}