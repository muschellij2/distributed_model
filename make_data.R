library(dplyr)
library(tibble)
library(readr)
set.seed(20190204)
n = 1000
n_sites = 3
indices = sample(1:n_sites, size = n, replace = TRUE)
true_beta = c(0.25, 1.25, -0.3)
df = tibble(
  ones = rep(1, n),
  x1 = rnorm(n),
  x2 = rnorm(n, mean = 0, sd = 2),
)
expb = exp(as.matrix(df) %*% true_beta )
df$prob_y = expb/(1 + expb)
df$y = rbinom(n, size = 1, prob = df$prob_y)
df$hospital = indices
sdf = split(df, df$hospital)
sdf = lapply(sdf, function(x) {
  x$hospital = NULL
  # random order this stuff
  x = x[ , sample(colnames(x))]
  x 
})
out_names = file.path("data", paste0("logistic_data_site", 1:n_sites, ".csv"))
mapply(function(x, y){
  readr::write_csv(x, y)
}, sdf, out_names)

df$hospital = NULL
formula = y ~ x1 + x2
full_model = glm(formula = formula, data = df, family = binomial())




my_binom = function(beta, df, formula) {
  y = model.frame(formula, data = df)[,1]
  X = model.matrix(formula, data = df)
  
  expb = exp(X %*% beta)
  p = expb / (1 + expb)
  p = c(p)
  gradient = colMeans(X * (y - p))
  new_beta = beta + gradient
  new_beta
}
beta = c(0.5, 0.5, 0.5)
new_beta = c(100, 100, 100)
for (i in 1:100) {
  new_beta = my_binom(beta, df, formula)
  if (max(abs(new_beta - beta)) <= 1e-5) {
    break
  }
  beta = new_beta
  print(i)
}
