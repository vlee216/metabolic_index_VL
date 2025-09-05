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
             taxa_id = g_i_i,
             minuslogpo2 = -log(all.dat$Pcrit),
             spc_in_PCgz = spc_in_PC_gz
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
    mu[i] <- Eo[taxa_id[i]]*invtemp[i] + n_pow[taxa_id[i]]*logW[i] - log(V[taxa_id[i]])
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

# Now for the method effect

# Fit model using method effect ####
## Setup TMB ####

data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i,
             minuslogpo2 = -log(all.dat$Pcrit),
             spc_in_PCgz = spc_in_PC_gz,
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
    mu[i] <- Eo[taxa_id[i]]*invtemp[i] + n_pow[taxa_id[i]]*logW[i] - log(V[taxa_id[i]] - method_mat[i,] %*% beta_method)
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
