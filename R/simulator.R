# library(OrvacRCT)


## a placeholder
test_mytrial <- function() {
  cat("nothing much here\n")
  test_trial1()
  test_trial2()
  test_trial3()
  test_trial4()
}



trial_dat <- function(){
  
  cfg <- readRDS("tools/cfg-example.RDS")
  cfg$people_per_month <- 50/3
  d <- get_trial_dat(cfg)
  d
}



