# library(OrvacRCT)
library(futile.logger)
library(optparse)



## a placeholder
test_mytrial <- function() {
  cat("nothing much here\n")
  test_trial1()
  test_trial2()
  test_trial3()
  test_trial4()
}






test_sim <- function(){
  
  cfg <- readRDS("R/cfg-example.RDS")
  d <- get_trial_dat(cfg)
  
}



get_cmdline_opts <- function(){
  
  option_list <- list(
    make_option(c("-f", "--cfgfile"), type = "character", default = "cfg1.yaml",
                help = "config file name", metavar = "character"),
    make_option(c("-o", "--use"), type = "logical", default = FALSE,
                help = "override config file with command line settings",
                metavar = "logical"),
    make_option(c("-i", "--idsim"), type = "character", default = FALSE,
                help = "label for current simulation", metavar = "character"),
    make_option(c("-n", "--nsims"), type = "integer", default = NULL,
                help = "number of simulations", metavar = "integer"),
    make_option(c("-s", "--seed"), type = "integer", default = NULL,
                help = "random seed", metavar = "integer"),
    make_option(c("-a", "--accrual"), type = "integer", default = NULL,
                help = "accrual rate i.e. people_per_interim_period",
                metavar = "integer"),
    make_option(c("-d", "--delay"), type = "double", default = NULL,
                help = "seroconversion information delay",
                metavar = "double"),
    make_option(c("-b", "--basesero"), type = "double", default = NULL,
                help = "baseline seroconversion prob", metavar = "double"),
    make_option(c("-p", "--trtprobsero"), type = "double", default = NULL,
                help = "trt seroconversion prob", metavar = "double"),
    make_option(c("-m", "--basemediantte"), type = "double", default = NULL,
                help = "baseline median time to med attendance (months)",
                metavar = "double"),
    make_option(c("-t", "--trtmedtte"), type = "double", default = NULL,
                help = "treatment arm median time to med attendance (months)",
                metavar = "double")
  )
  
  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);
  opt
  
}


get_cfg <- function(cfgfile = "cfg1.yaml"){
  
  opt <- get_cmdline_opts()
  
  cat(opt)
  
  tt <- tryCatch(configtmp <- read.config(file = file.path(getwd(), "tools", cfgfile)),
                 error=function(e) e,
                 warning=function(w) w)
  ifelse(is(tt,"warning") | is(tt,"error"),"Configuration Warning/Error.
         Please ensure configuration file has terminating empty line.",
         "Configuration File Loaded OK")
  
  l <- list()
  
  #l$dnames <- dnames

  l$n_sims <- tt$n_sims
  l$n_start <- tt$n_start
  l$n_stop <- tt$n_stop
  l$accrual_per_qtr <- tt$accrual_per_qtr
  l$sero_info_delay <- tt$sero_info_delay
  l$n_max_sero <- tt$n_max_sero
  l$n_start_clin <- tt$n_start_clin
  l$age_months_lwr <- tt$age_months_lwr
  l$age_months_upr <- tt$age_months_upr
  l$max_age_fu_months <- tt$max_age_fu_months
  l$baseline_prob_sero <- tt$baseline_prob_sero
  l$trtprobsero <- tt$trtprobsero
  l$deltaserot3 <- compute_sero_delta(l$baseline_prob_sero, l$trtprobsero)
  
  
  # colourblind palette
  l$cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  return(l)
  
  
}




compute_sero_delta <- function(p_ctl, p_trt){
  
  delta <- (p_trt - p_ctl)/ (1 - p_ctl)
  delta
  
}
