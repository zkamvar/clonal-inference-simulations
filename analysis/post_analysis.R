library('zksimanalysis')
library('magrittr')
library('viridis')
library('dplyr')
library('tidyr')
library('purrr')
library('ggplot2')

datfiles <- dir("rda_files/", full.names = TRUE) %>% grep("feather.rda$", ., value = TRUE)

for (i in datfiles){
  cat(i, "\n")
  load(i)
}

ex_run  <- "twenty_loci([0-9]+?)/"
ex_seed <- "seed_([0-9]+?)_"
ex_sex  <- "sex_([0-9.]+?)_"
ex_gen  <- "gen_([0-9]+?)_"
ex_rep  <- "rep_([0-9]+?).pop"
ex_samp <- "_sam_([0-9]+?)$"

datnames <- grep("^X.+?_twenty.+?feather$", ls(), value = TRUE)

datalist <- datnames %>%
  lapply(get) %>%
  setNames(datnames) %>%
  bind_rows(.id = "source") %>%
  select(Ia:pop, source) %>%
  extract(pop, c("run", "seed", "sexrate", "gen", "rep", "sample"),
          paste0(ex_run, ex_seed, ex_sex, ex_gen, ex_rep, ex_samp),
          remove = FALSE) %>%
  mutate(sample = factor(sample, unique(sample)))

vals <- datalist %>%
  mutate(mean.rd = vapply(samples.rd, mean, numeric(1), na.rm = TRUE)) %>%
  mutate(sd.rd = vapply(samples.rd, sd, numeric(1), na.rm = TRUE)) %>%
  mutate(var.rd = vapply(samples.rd, var, numeric(1), na.rm = TRUE)) %>%
  mutate(min.rd = vapply(samples.rd, min, numeric(1), na.rm = TRUE)) %>%
  mutate(max.rd = vapply(samples.rd, max, numeric(1), na.rm = TRUE)) %>%
  select(-starts_with("samples"))

ggplot(vals, aes(x = sexrate, y = rbarD, fill = sample)) +
  geom_boxplot()

ggplot(vals, aes(x = sexrate, y = p.rD, fill = sample)) +
  geom_boxplot() +
  scale_y_log10(breaks = c(0.01, 0.1, 1),
                minor_breaks = c((1:9)*0.001, (2:9)*0.01, (2:9)*0.1))

ggplot(vals, aes(x = sexrate, y = rbarD, color = log(p.rD))) +
  geom_point(alpha = 0.25, position = "jitter") +
  geom_boxplot(color = "black", notch = TRUE, width = 0.5, alpha = 0.25) +
  scale_color_viridis(option = "viridis", breaks = log(c(0.005, 0.01, 0.025, 0.05, 0.1)), labels = c(0.005, 0.01, 0.025, 0.05, 0.1)) +
  facet_wrap(~sample, nrow = 1)

p1 <- vals %>%
  filter(p.rD == 1, !is.na(rbarD))



