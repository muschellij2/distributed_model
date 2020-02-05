for (i in 1:100) {
  print(i)
  source("secondary.R")
  source("master.R")
}
