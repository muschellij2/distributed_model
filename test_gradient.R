library(readr)
df = read_csv("data/logistic_data_site1.csv")
shuffle_rows = TRUE
beta = rep(0, 3)
family = binomial()
formula = y ~ x1 + x2
gradient_value = function(beta = NULL, df, formula, 
                          family = binomial(), 
                          iteration_number = 0,
                          shuffle_rows = TRUE) {
  if (shuffle_rows) {
    df = df[sample(nrow(df)), ]
  }
  y = model.frame(formula, data = df)[,1]
  X = model.matrix(formula, data = df)
  cc = complete.cases(X)
  X = X[cc,]
  y = y[cc]
  
  if (is.null(beta)) {
    beta = rep(0, ncol(X))
  }
  linkinv = family$linkinv
  variance <- family$variance
  mu.eta <- family$mu.eta
  dev.resids <- family$dev.resids
  eta = drop(X %*% beta)
  mu = linkinv(eta)
  varmu <- variance(mu)
  mu.eta.val <- mu.eta(eta)
  weights = rep(1, length = length(eta))
  good <- (weights > 0) & (mu.eta.val != 0)  
  offset = 0
  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
  w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
  xwtx = t(X) %*% diag(w) %*% X
  # dev <- sum(dev.resids(y, mu, weights))
  # w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
  # residuals <- (y - mu)/mu.eta(eta)  
  
  
  p = linkinv(X %*% beta)
  # expb = exp(X %*% beta)
  # p = expb / (1 + expb)
  p = c(p)
  gradient = colSums(X * (y - p))
  stopifnot(length(gradient) == length(beta))
  result = list(
    gradient = gradient,
    sample_size = nrow(X),
    iteration_number = iteration_number
  )
  return(result)
}
