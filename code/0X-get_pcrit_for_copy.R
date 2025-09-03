#options(echo = FALSE)
rm(list = ls())
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(worrms)
library(purrr)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Setup Code ####
## load functions ####
source("code/helper/fit_model_funs.R")

# load (or create) fitted model
# Future testing note -- if you have changed anything that needs updated, set below to = T instead of = F
model.fit <- fit_model_augmented_taxa(fitnew = F)
model.fit$rep

# Do this for several types of salmonidae as a test
species.names <- c("Oncorhynchus mykiss", "Oncorhynchus nerka", "Salvelinus albus")

# Set method names according to the leveling in 02-fit_tmb_model
method.names <- c("oxyconform", "death", "smr")

# create a dataframe of all combinations of the above
prediction_info <- tidyr::expand_grid(taxa.name = species.names, method = method.names)


# calcualte pcrit for each combination above and put in data frame
# Now that I have multiple lines for beta_method, this throws an error. Will need editing in fit_model_funs.R
result_df <- pmap_dfr(prediction_info, 
                      function(taxa.name, method) {
                        result <- estimate_taxa(taxa.name = taxa.name,
                                                w = 25,
                                                temperature = 20,
                                                method = method,
                                                rep = model.fit$rep,
                                                ParentChild_gz = model.fit$ParentChild_gz,
                                                ps = 0)
                        tibble(species = taxa.name,
                               method = method,
                               log_pcrit = as.numeric(result$log_pcrit["logpcrit"]),
                               se_pcrit = as.numeric(result$log_pcrit["se"] )
                        )
                        }
                      )


# rename species for genus and family - level
result_df$species[result_df$species == species.names[2] ] <- "Oncor sp."
result_df$species[result_df$species == species.names[3] ] <- "Salmonidae"

# adjust to allow for plotting smr and routine jittered vertically
result_df <- result_df %>%
  mutate(
    species_factor = factor(species, levels = rev(unique(species))),  # so species appear in the same order as the plot
    ypos = as.numeric(species_factor) + ifelse(method == "smr", -0.1, 0.1)
  )

result_df

# Current issue -- the SE values for Cottidae are smaller than should be expected
# They should be the largest SE values for each of the methods, not the smallest
# What's going on? did the dataset acquire a Astrocottus leprops value at some point?

# colors to use:
colors = c("#0072B2", "#D55E00", "#ADC905")
# setup theme
theme_set(theme_bw(base_size = 18))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.text = element_text(color = "black"),
             legend.text = element_text(size = 16)
             )
# make plot
plot_results <- ggplot(data = result_df, aes(y = ypos, x = log_pcrit, col = method)) + 
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = log_pcrit - se_pcrit,
                     xmax = log_pcrit + se_pcrit),
                 linewidth = 1, height = 0.1,
                 show.legend = FALSE) +
  scale_y_continuous(
    breaks = 1:length(unique(result_df$species_factor)),
    labels = levels(result_df$species_factor)
  ) +
  scale_x_continuous(
    breaks = log(seq(2, 10, by = 1)),  
    labels = seq(2, 10, by = 1),
    minor_breaks = log(seq(2,10, by = 0.1)),
    name = bquote(p[crit]) ) +
  ylab("") +
  scale_color_manual(values = c("death" = colors[1], "smr" = colors[2], "oxyconform" = colors[3]),
                     name = NULL,
                     breaks = c("death", "smr", "oxyconform"),
                     labels = c("Death", "SMR", "Oxyconform")) +
  theme(legend.position = "top")

print(plot_results)
ggsave(file = "figures/salmonidae.png",
       plot = plot_results,
       device = "png",
       units = "px",
       scale = 5,
       width = 600,
       height = 300)

# As expected -- unfortunately -- issue persists for Salmonidae test family
# The SE issue is not being caused by the choice of family
# This is coming from the some kind of code or model intrinsic issue
