#' ---
#' title: "Assessing bias due to reduced sample size"
#' output:
#'   html_document: default
#'   html_notebook: default
#' ---
#'
#' ## Introduction
#'
#' We've shown that clone-correction can reduce the power of $\bar{r}_d$ to detect
#' clonal reproduction, but it's also clear that a reduced sample size contributes
#' to this as well. To disentangle this, we assess the $\bar{r}_d$ from 1000
#' replicates of *n* samples where *n* is the size of the data set after
#' clone-correction. The data were processed via the `analyze_jackknife_ia.R`
#' script. These were saved in the `ma_jack_rda_files/` directory.
#'
#'
#' ## Setup
#'
#'
#' ### Packages
#'
## ---- required_packages--------------------------------------------------
library('zksimanalysis')

#'
#' ### Loading Data
#'
## ------------------------------------------------------------------------
res        <- load("data/processed_results.rda") # load original values of rd
jack_files <- list.files("data/jack_rda_files/", full.names = TRUE) # load jackknife simulations
jacklist   <- setNames(jack_files, basename(jack_files))

for (i in jack_files){
  jacklist[basename(i)] <- load(i)
}
jackdat <- jacklist %>%
  map_df(get, .id = "source") %>%
  mutate(mutation_rate = ifelse(grepl("mutant", source), "even", "uneven")) %>%
  mutate(source = gsub("(^.+?feather).*$", "\\1.DATA.rda", source)) %>%
  select(matches("jack"), everything())
rm(list = jacklist)
gc()

# taking stock of how much data we have
jackdat %>% filter(!is.na(pop))
vals

#'
#' Note that above, I'm filtering `jackdat` for non-missing populations. I should
#' explain this a bit. When `jack.ia` was put into *poppr*, the requirements for
#' running that **n** was greater than 2 and smaller than the total popualtion size.
#' If the number of MLGs did not fit this requirement, `tidy_jackia` would  return
#' all `NA`, including that for population.
#'
#' We can't ask any questions regarding bias with missing data, so we are filtering
#' it out here.
#'
## ------------------------------------------------------------------------
# combining original values and jackknife simulations
vals <- get(res) %>%
  inner_join(filter(jackdat, !is.na(pop)), # filtering out missing data
             by = c("source", "mutation_rate", "pop")) %>%
  select(rbarD, jack.rd, jack.p.rD, mutation_rate, everything())
vals
even_mutation   <- vals %>% filter(mutation_rate == "even")
uneven_mutation <- vals %>% filter(mutation_rate == "uneven")

#'
#' ## Data Exploring!
#'
#' First, let's see what we have after the filtering:
#'
## ------------------------------------------------------------------------
vals %>%
  group_by(mutation_rate, sexrate, sample) %>%
  summarize(n = n()) %>%
  spread(sample, n) %>%
  knitr::kable()

#' Now, let's visualize the "p-value" for each jack knife analysis. This p-value
#' represents the fraction of jack-knife simulations less than or equal to the
#' clone-corrected value.
#'
## ------------------------------------------------------------------------
vals %>%
  filter(p.rD > 0.01, sexrate < 0.001) %>%
  ggplot(aes(x = jack.p.rD, fill = mutation_rate)) +
  geom_histogram(alpha = 0.5) +
  facet_grid(sexrate~sample, scale = "free_y")

#'
