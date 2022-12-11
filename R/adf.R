"adf" <- function(x,k = 0, int = TRUE, trend = FALSE){
  require(dynlm)
  dx <- diff(x)
  formula <- paste("dx ~ L(x)")
  if(k > 0)
    formula <- paste(formula," + L(dx,1:k)")
  if(trend){
    s <- time(x)
    t <- ts(s - s[1],start = s[1],freq = frequency(x))
    formula <- paste(formula," + t")
  }
  if(!int) formula <- paste(formula," - 1")
  summary(dynlm(as.formula(formula)))
}
