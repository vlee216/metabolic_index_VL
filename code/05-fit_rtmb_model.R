# Code fit hierarchical model to Arrenhius Equation
rm(list = ls())
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(mvtnorm)
library(tmbstan)
library(shinystan)
library(egg)
library(cowplot)
library(rlang)

# Setup Code ####


## load functions ####
source("code/helper/fit_model_funs.R")
source("code/helper/calculate_EDF_fn.R")

# ggplot setup ####
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.text = element_text(color = "black"),
             axis.text.y = element_text(size = 16)
)
# ploting and summary functions used below ####
make_group_plot <- function(group_to_plot,sum_est, saveplot = T) {
  group_sym <- enquo(group_to_plot)
  dataset_name <- paste0(as_name(group_sym), "Est")
  
  Est.2.plot <- sum_est[[dataset_name]] %>%
    filter(NoSpecies  >= 2) %>%
    mutate(group_num = paste0(!!group_sym, " (",  NoSpecies , ")"))
  
  Vplot <- plotest(Est.2.plot, logV,  logVse, group_num)
  Vplot <- Vplot + xlab("log(V)") + xlim(x_v_lims) + ylab("")
  Eoplot <- plotest(Est.2.plot, Eo, Eose, group_num)
  Eoplot <- Eoplot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())  + 
    xlab(expression("E"[o])) + xlim(x_eo_lims) 
  
  nplot <- plotest(Est.2.plot, n, nse, group_num)
  nplot <- nplot + xlim(x_n_lims) + theme(axis.text.y = element_blank(), 
                                          axis.title.y = element_blank()
  )
  
  group_plot <- ggarrange(Vplot, 
                          nplot,
                          Eoplot,
                          nrow = 1)
  print(group_plot)
  if (saveplot) {
    saveplotname <- paste0("figures/", as_name(group_sym), "_plot_mle.png")
    ggsave(filename = saveplotname,
           plot = group_plot,
           device = "png",
           units = "px",
           scale = 3,
           width = 1029,
           height = 629)
  }
  return(group_plot)
}
# Make a dataframe for plotting group-level estimates #
make_df_plot <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  Est <- tibble("{groupname}" := GroupNames,
                logV = beta_mle[group_index,1],
                logVse =  beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                Eose = beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nse =  beta_se[group_index,2]
  )
  return(Est)
}
# Function for plotting estimates by specified group ####
plotest <- function(dataest, trait, xse, groupname) {
  trait <- enquo(trait)
  groupname <- enquo(groupname)
  xse <- enquo(xse)
  
  groupplot <- ggplot(data = dataest, aes(x = !!trait, y = !!groupname)) +
    geom_point(size = 2) + 
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(y = !!groupname,
                      xmin = !!trait - !!xse,
                      xmax = !!trait + !!xse)
    )
  
  return(groupplot)
}
# Function to make a species-level dataframe of estimates ####
make_species_df <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  
  Est <- tibble("{groupname}" := GroupNames,
                logV = beta_mle[group_index,1],
                logVSE = beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                EoSE = beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nSE = beta_se[group_index,2],
                
  )
  return(Est)
}
summarize_estimates <- function (beta_mle, beta_se, ParentChild_gz, taxa.list){
  
  if (taxa.list[1] == "Family") baselevel = 0
  if (taxa.list[1] == "Order") baselevel =1
  if (taxa.list[1] == "Class") baselevel =2
  if (taxa.list[1] == "Phylum") baselevel =3
  
  ## Phylum ####
  if (taxa.list[1] == "Phylum") {
    PhylumEst <- make_df_plot(level = 1,
                              beta_mle,
                              beta_se,
                              ParentChild_gz,
                              groups = taxa.list)
    PhylumEst <- merge(PhylumEst, phylum_summary)
  }
  ## Class ####
  if (taxa.list[1] %in% c("Phylum", "Class")) {
    ClassEst <- make_df_plot(level = 2,
                             beta_mle,
                             beta_se,
                             ParentChild_gz,
                             groups = taxa.list)
    ClassEst <- merge(ClassEst, class_summary)
  }
  if (taxa.list[1] %in% c("Order", "Class", "Phylum") ) {
    ## By Order #####
    OrderEst <- make_df_plot(level = baselevel, 
                             beta_mle,
                             beta_se,
                             ParentChild_gz,
                             groups = taxa.list)
    
    OrderEst <- merge(OrderEst, order_summary)
  }
  if (taxa.list[1] %in% c("Family","Order", "Class", "Phylum")) {
    ## By Family ####
    FamilyEst <- make_df_plot(level = baselevel + 1, 
                              beta_mle,
                              beta_se,
                              ParentChild_gz,
                              groups = taxa.list)
    FamilyEst <- merge(FamilyEst, family_summary)
  }
  
  
  SpeciesEst <- make_species_df(level = baselevel + 3, 
                                beta_mle,
                                beta_se,
                                ParentChild_gz,
                                groups = taxa.list)
  
  
  output <- list(FamilyEst = FamilyEst,
                 SpeciesEst = SpeciesEst)
  if (taxa.list[1] == "Phylum") output$PhylumEst = PhylumEst
  if (taxa.list[1] %in% c("Phylum", "Class") ) output$ClassEst = ClassEst
  if (taxa.list[1] %in% c("Phylum", "Class", "Order") ) output$OrderEst = OrderEst
  return(output)
}

# load data ####
all.dat <- load_data()
all.dat <- filter_data(all.dat)

# Transform EstMethod_Metric == "loe" to EstMethod_Metric == "lethal (VL, 8/27)
for (i in 1:nrow(all.dat)){
  if (all.dat[i,]$EstMethod_Metric == "loe") {
    all.dat[i,]$EstMethod_Metric <- "death"
  }
}

# Filter to remove EstMethod_Metric that is not "oxyconform", "smr", or "death"
all.dat <- all.dat %>%
  dplyr::filter(EstMethod_Metric == "oxyconform" | EstMethod_Metric == "smr" | EstMethod_Metric == "death")

# This should produce TRUE
length(unique(all.dat$EstMethod_Metric)) == 3

# change so that the default level is oxyconform 
all.dat$EstMethod_Metric <- as.factor(all.dat$EstMethod_Metric)
all.dat$EstMethod_Metric <- factor(all.dat$EstMethod_Metric, levels = c("oxyconform", "smr", "death"))

method_mat <- model.matrix(~ EstMethod_Metric , all.dat)
n_methods <- ncol(method_mat) -1

 # Summarise data by taxa ####
phylum_summary <- all.dat %>%
  group_by(Phylum) %>%
  summarise(NoClass = n_distinct(Class), NoOrder = n_distinct(Order), NoFamily = n_distinct(Family),  NoGenus = n_distinct(Genus), NoSpecies = n_distinct(Species))

class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genus)), NoSpecies= length(unique(Species)))

genus_summary <- all.dat %>%
  group_by(Family, Order, Genus) %>%
  summarise(NoSpecies= length(unique(Species)))

# Setup data for TMB ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Phylum", "Class","Order", "Family", "Genus", "Species")
### Create new ParentChild matrix for reduced taxonomic structure ####
taxa.info <- make_taxa_tree(all.dat, taxa.list)
ParentChild_gz <- taxa.info$ParentChild_gz
PC_gz <- taxa.info$PC_gz
g_i <- taxa.info$g_i
g_i_i <- taxa.info$g_i_i
n_k <- taxa.info$n_k
n_j <- taxa.info$n_j
n_g <- taxa.info$n_g
n_i <- taxa.info$n_i
spc_in_PC_gz <- taxa.info$spc_in_PC_gz

# Fit model without method effect ####
## Setup TMB ####
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = -log(all.dat$Pcrit),
             spc_in_PCgz = spc_in_PC_gz -1
)

parameters = list(alpha_j = c(0, 0, 0),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  logsigma = 0
)

## space to work on turning this into RTMB code

# load the package
library(RTMB)

# build the cov_matrix() function
cov_matrix <- function(L_val, min_var, n_rows){
  L_rc <- matrix(0, nrow = n_rows, ncol = n_rows)
  Cov_rr <- matrix(0, nrow = n_rows, ncol = n_rows)
  Return_rr <- matrix(0, nrow = n_rows, ncol = n_rows)
  Count <- 1
  for (r in 1:n_rows){
    for (c in 1:n_rows){
      if (r>=c){
        L_rc[r,c] <- L_val[Count]
        Count <- Count + 1
      }
      else{
        L_rc[r,c] <- 0.0
      }
    }
  }
  diag(Cov_rr) <- min_var
  Cov_rr <- Cov_rr + L_rc %*% t(L_rc)
  Return_rr <- Cov_rr
  Return_rr
}

# set up the RTMB model to mirror the previous C++ code

f_base <- function(parameters){
  # unpack the data and parameters
  getAll(parameters, data)
  
  # establish the internal values taken from the provided data/parameters
  sigma <- exp(logsigma)
  n_j <- ncol(beta_gj)
  n_i <- length(spc_in_PCgz)
  n_g <- nrow(PC_gz)
  n_d <- length(minuslogpo2)
  n_l <- length(log_lambda)
    
  # build needed vectors and matrices
  Parent_j <- vector(mode = "numeric", length = n_j)
  Prediction_j <- vector(mode = "numeric", length = n_j)
  Deviation_j <- vector(mode = "numeric", length = n_j)
  tmpCov_jj <- matrix(0, nrow = n_j, ncol = n_j)
  V <- vector(mode = "numeric", length = n_i)
  n_pow <- vector(mode = "numeric", length = n_i)
  Eo <- vector(mode = "numeric", length = n_i)
  spc_ij <- matrix(0, nrow = n_i, ncol = n_j)
  lambda <- vector(mode = "numeric", length = n_l)
  
  # build the spc_ij matrix
  for (i in 1:n_i){
    for (j in 1:n_j){
      if (j == 1){
        spc_ij[i,j] <- exp(beta_gj[spc_in_PCgz[i],j])
      }
      if (j > 1){
        spc_ij[i,j] <- beta_gj[spc_in_PCgz[i],j]
      }
    }
  }
  
  # make vector of lambdas
  lambda <- exp(log_lambda)
  
  # extract more needed values
  V <- spc_ij[,1]
  n_pow <- spc_ij[,2]
  Eo <- spc_ij[,3]
  mu <- vector(mode = "numeric", length = n_d)
  
  # initialize the jnll and set to zero
  jnll_comp <- c(rep(0,3))
  
  # set priors on log_lambda
  df <- 3
  tsigma <- 5
  nlambda <- length(log_lambda)
  jnll_comp[3] <- -sum(dt(log_lambda, df, log = TRUE)) - log(tsigma)*nlambda
    
  # build covariance matrix using the cov_matrix() established function
  min_var <- 0.001
  Cov_jj <- matrix(0, nrow = n_j, ncol = n_j)
  Cov_jj <- cov_matrix(L_z, min_var, n_j)
  
  # random effects
  covmult <- NA
  for (g in 1:n_g){
    Child_num <- PC_gz[g,2]
    Parent_row <- PC_gz[g,1] + 1
    lambda_num <- Child_num
    for (j in 1:n_j){
      if(PC_gz[g,2] == 0) {Parent_j[j] <- alpha_j[j]}
      if(PC_gz[g,2] >= 1) {Parent_j[j] <- beta_gj[Parent_row, j]}
      Prediction_j[j] <- Parent_j[j]
    }
    for (j in 1:n_j){
      Deviation_j[j] <- beta_gj[g,j] - Prediction_j[j]
    }
    if( Child_num == 0) {covmult <- 1}
    if( Child_num >= 1) {covmult <- lambda[lambda_num]}
    tmpCov_jj <- Cov_jj * covmult
    jnll_comp[1] <- jnll_comp[1] - dmvnorm(Deviation_j, Sigma = tmpCov_jj, log = TRUE)
  }
  
  # probability of the data
  for (i in 1:n_d){
    mu[i] <- Eo[taxa_id[i] + 1]*invtemp[i] + n_pow[taxa_id[i] + 1]*logW[i] - log(V[taxa_id[i] + 1])
  }
  
  # get nll of the data
  jnll_comp[2] <- -sum(dnorm(minuslogpo2, mu, sigma, log = TRUE))
  
  # provide total jnll
  jnll <- sum(jnll_comp)
  nll_data <- jnll_comp[2]
  
  # adreport / report
  ADREPORT(spc_ij)
  REPORT(Cov_jj)
  REPORT(mu)
  REPORT(nll_data)

  # print jnll
  jnll
}


# complete model fit using RTMB
obj_base <- MakeADFun(f_base, parameters, random = "beta_gj")
opt_base <- nlminb(obj_base$par, obj_base$fn, obj_base$gr)

sdr <- sdreport(obj_base)

sdr
# not sure yet if this is working right -- the code runs, but is 
# this giving me the model I'm supposed to get?
# will interrogate this over the next workweek

## end space to work on turning this into RTMB code


Random <- c("beta_gj")
model <- "hierarchical_mi_base"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj_nomethod <-
  TMB::MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )

opt_nomethod <- nlminb(obj_nomethod$par, obj_nomethod$fn, obj_nomethod$gr)
rep_nomethod = TMB::sdreport( obj_nomethod,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)
EDF <- calculate_EDF( obj= obj_nomethod,
                      opt = opt_nomethod,
                      nonvariance_fixed_effects = c( "alpha_j", "L_z", "log_lambda"),
                      prediction_name = "mu",
                      data_name = "minuslogpo2",
                      delta = 0.01,
                      show_progress = F,
                      refit = "random"
)

nll_data <- obj_nomethod$report()[["nll_data"]]
cAIC_nomethod <- 2 * nll_data + 2 * EDF


# Now for the method effect

# Fit model using method effect ####
## Setup TMB ####

data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = -log(all.dat$Pcrit),
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1]
)

parameters = list(alpha_j = c(0, 0, 0),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = rep(0, ncol(method_mat) -1),
                  logsigma = 0
)


## space for RTMB code (method)



f <- function(parameters){
  # unpack the data and parameters
  getAll(parameters, data)
  
  # establish the internal values taken from the provided data/parameters
  sigma <- exp(logsigma)
  n_j <- ncol(beta_gj)
  n_i <- length(spc_in_PCgz)
  n_g <- nrow(PC_gz)
  n_d <- length(minuslogpo2)
  n_l <- length(log_lambda)
  
  # build needed vectors and matrices
  Parent_j <- vector(mode = "numeric", length = n_j)
  Prediction_j <- vector(mode = "numeric", length = n_j)
  Deviation_j <- vector(mode = "numeric", length = n_j)
  tmpCov_jj <- matrix(0, nrow = n_j, ncol = n_j)
  V <- vector(mode = "numeric", length = n_i)
  n_pow <- vector(mode = "numeric", length = n_i)
  Eo <- vector(mode = "numeric", length = n_i)
  spc_ij <- matrix(0, nrow = n_i, ncol = n_j)
  lambda <- vector(mode = "numeric", length = n_l)
  
  # build the spc_ij matrix
  for (i in 1:n_i){
    for (j in 1:n_j){
      if (j == 1){
        spc_ij[i,j] <- exp(beta_gj[spc_in_PCgz[i],j])
      }
      if (j > 1){
        spc_ij[i,j] <- beta_gj[spc_in_PCgz[i],j]
      }
    }
  }
  
  # make vector of lambdas
  lambda <- exp(log_lambda)
  
  # extract more needed values
  V <- spc_ij[,1]
  n_pow <- spc_ij[,2]
  Eo <- spc_ij[,3]
  mu <- vector(mode = "numeric", length = n_d)
  
  # initialize the jnll and set to zero
  jnll_comp <- c(rep(0,3))
  
  # set priors on log_lambda
  df <- 3
  tsigma <- 5
  nlambda <- length(log_lambda)
  jnll_comp[3] <- -sum(dt(log_lambda, df, log = TRUE)) - log(tsigma)*nlambda
  
  # build covariance matrix using the cov_matrix() established function
  min_var <- 0.001
  Cov_jj <- matrix(0, nrow = n_j, ncol = n_j)
  Cov_jj <- cov_matrix(L_z, min_var, n_j)
  
  # random effects
  covmult <- NA
  for (g in 1:n_g){
    Child_num <- PC_gz[g,2]
    Parent_row <- PC_gz[g,1] + 1
    lambda_num <- Child_num
    for (j in 1:n_j){
      if(PC_gz[g,2] == 0) {Parent_j[j] <- alpha_j[j]}
      if(PC_gz[g,2] >= 1) {Parent_j[j] <- beta_gj[Parent_row, j]}
      Prediction_j[j] <- Parent_j[j]
    }
    for (j in 1:n_j){
      Deviation_j[j] <- beta_gj[g,j] - Prediction_j[j]
    }
    if( Child_num == 0) {covmult <- 1}
    if( Child_num >= 1) {covmult <- lambda[lambda_num]}
    tmpCov_jj <- Cov_jj * covmult
    jnll_comp[1] <- jnll_comp[1] - dmvnorm(Deviation_j, Sigma = tmpCov_jj, log = TRUE)
  }
  
  # probability of the data
  for (i in 1:n_d){
    mu[i] <- Eo[taxa_id[i] + 1]*invtemp[i] + n_pow[taxa_id[i] + 1]*logW[i] - log(V[taxa_id[i] + 1] - method_mat[i,] %*% beta_method)
  }
  
  # get nll of the data
  jnll_comp[2] <- -sum(dnorm(minuslogpo2, mu, sigma, log = TRUE))
  
  # provide total jnll
  jnll <- sum(jnll_comp)
  nll_data <- jnll_comp[2]
  
  # adreport / report
  ADREPORT(spc_ij)
  ADREPORT(lambda)
  REPORT(Cov_jj)
  REPORT(mu)
  REPORT(nll_data)
  
  # print jnll
  jnll
}


# complete model fit using RTMB
obj <- MakeADFun(f, parameters, random = "beta_gj")
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)

sdr
# the method RTMB code is not currently functional -- these estimates are completely off.
# wrong sign
# they need to be negative not positive

## end space for RTMB code (method)



Random <- c("beta_gj")
model <- "hierarchical_mi"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj <-
  TMB::MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = TMB::sdreport( obj,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)

## calculate cAIC
EDF <- calculate_EDF( obj= obj,
                      opt = opt,
                      nonvariance_fixed_effects = c("beta_method", "alpha_j", "L_z", "log_lambda"),
                      prediction_name = "mu",
                      data_name = "minuslogpo2",
                      delta = 0.01,
                      show_progress = F,
                      refit = "random"
)

nll_data <- obj$report()[["nll_data"]]
cAIC_method <- 2 * nll_data + 2 * EDF
cAIC_list <- c(cAIC_nomethod, cAIC_method)
print(cAIC_list - min(cAIC_list))


re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)
beta_method <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),1], nrow = n_methods)
trans <- summary(rep, "report")
lambda <- trans[grep(rownames(trans), pattern = "lambda"), ]


# Summarize Estimatess ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

# Plot Estimates ####
# set axis limits
x_v_lims <- c(0.75, 2.2)
x_eo_lims <- c(-0.125, 0.6)
x_n_lims <- c(-0.2, 0.075)

plot_est <- T
if (plot_est) {
  if (taxa.list[1] == "Phylum") make_group_plot(Phylum, sum_est, saveplot = T)

  if (taxa.list[1] %in% c("Phylum", "Class") ) make_group_plot(Class, sum_est, saveplot = T)
  make_group_plot(Order, sum_est, saveplot = T)
  make_group_plot(Family, sum_est,  saveplot = T)
}
# Compare fits to individual species ####
SpeciesEst <- make_species_df(level = which(taxa.list == "Species"), 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

SpeciesEst$V = exp(SpeciesEst$logV)
saveRDS(SpeciesEst, file = "analysis/hierarchical_species_estimates.RDS")

p_diagnostic <- plot_diagnostics(model= "method", 
                                 Pcrit = all.dat$Pcrit, 
                                 inv.temp = all.dat$inv.temp, 
                                 W = all.dat$W, 
                                 SpeciesEst = SpeciesEst,
                                 beta_method = beta_method,
                                 method_mat = method_mat)
ggsave(filename= "figures/diagnostic.png",
       plot = p_diagnostic,
       units = "px",
       scale = 3,
       width = 1029,
       height = 1029)


# Get species-level trait variance and make table
sampled_species_standard_deviation <- apply(X = SpeciesEst[,c("logV","n", "Eo")], MAR = 2, FUN = sd)



## a plot of all species pcrit, vs. temperature and body mass
temp_2_plot <- seq(0, 35, by = 5)
inv_temp_ticks <- (1/ kb) * (1 / kelvin(temp_2_plot) - 1 / kelvin(tref) )

alldata_plot_temperature <- ggplot(data = all.dat, aes(x = inv.temp, y = log(Pcrit), col = as.factor(Phylum) )) + 
  scale_colour_viridis_d(option = "turbo") +
  geom_point(size = 2) +
  scale_x_continuous(breaks = inv_temp_ticks, labels = temp_2_plot,
                     name = "Temperature °C") +
  ylab(bquote( log( p["crit"] ) ) ) +
  labs(col = "Phylum")

alldata_plot_w<- ggplot(data = all.dat, aes(x =log(W), y = log(Pcrit), col = as.factor(Phylum) )) + 
  scale_colour_viridis_d(option = "turbo") +
  geom_point(size = 2) +
  xlab ("log(W)" ) + 
  ylab(bquote( log( p["crit"] ) ) ) +
  labs(col = "Phylum")

alldata_plot <- 
  cowplot::plot_grid(alldata_plot_temperature + theme(legend.position="none"),
             alldata_plot_w + theme(legend.position="none"),
             nrow = 1,
             align = "v")
  
legend <- get_legend(
    # create some space to the left of the legend
    alldata_plot_temperature + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
alldata_plot <- plot_grid(alldata_plot, legend, rel_widths = c(3, 0.85) )

ggsave(file = "figures/alldata_plot.png",
       plot = alldata_plot,
       units = "px",
       scale = 2,
       height = 500,
       width = 1200,
       bg = "white")

# print out parameter values reported in manuscript

## levels of variance - the log_lambdas
print(exp(fixef[rownames(fixef) == "log_lambda", "Estimate"]))

## print out alpha_j - the mean trait values
fixef[rownames(fixef) == "alpha_j", ]
### First element is log(V).  Convert to V and get SE using delta method
c(exp(fixef[rownames(fixef) == "alpha_j", "Estimate"][1])  , se =  exp(fixef[rownames(fixef) == "alpha_j", "Estimate"][1]) * fixef[rownames(fixef) == "alpha_j", "Std. Error"][1] )

# for each species, document how many unique body sizes and unique temperatures

number_bodysize_temperature <- all.dat %>%
  group_by(Species) %>%
  summarise(n_sizes = length(unique(W)), n_temps = length(unique(Temp)))
n_species_2_temperature <- length(which(number_bodysize_temperature$n_temps > 2))
n_species_2_sizes <- length(which(number_bodysize_temperature$n_sizes > 2))
print(c(n_species_2_sizes, n_species_2_temperature))

