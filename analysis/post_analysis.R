library('zksimanalysis')
library('magrittr')
library('viridis')
library('dplyr')
library('tidyr')
library('purrr')
library('ggplot2')

datfiles <- dir("rda_files/", full.names = TRUE) %>% grep("feather.rda$", ., value = TRUE)
divfiles <- dir("diversity_rda_files/", full.names = TRUE)
locfiles <- dir("locus_rda_files/", full.names = TRUE)


for (i in datfiles){
  cat(i, "\n")
  load(i)
}

for (i in divfiles){
  cat(i, "\n")
  load(i)
}

for (i in locfiles){
  cat(i, "\n")
  load(i)
}


ex_run  <- "twenty_loci([0-9]+?)/"
ex_seed <- "seed_([0-9]+?)_"
ex_sex  <- "sex_([0-9.]+?)_"
ex_gen  <- "gen_([0-9]+?)_"
ex_rep  <- "rep_([0-9]+?).pop"
ex_samp <- "_sam_([0-9]+?)$"

datnames <- grep("^X.+?twenty.+?feather$", ls(), value = TRUE)
divnames <- grep("^X.+?twenty.+?feather.divtable.rda$", ls(), value = TRUE)
locnames <- grep("^X.+?twenty.+?feather.locustable.rda$", ls(), value = TRUE)

datalist <- datnames %>%
  lapply(get) %>%
  setNames(datnames) %>%
  bind_rows(.id = "source") %>%
  select(Ia:rbarDcc, source) %>%
  extract(pop, c("run", "seed", "sexrate", "gen", "rep", "sample"),
          paste0(ex_run, ex_seed, ex_sex, ex_gen, ex_rep, ex_samp),
          remove = FALSE) %>%
  mutate(sample = factor(sample, unique(sample)))

divlist <- divnames %>%
  lapply(get) %>%
  setNames(divnames) %>%
  bind_rows()

loclist <- locnames %>%
  lapply(get) %>%
  setNames(locnames) %>%
  bind_rows()


CF <- function(NMLG, sample){
  sample <- as.integer(as.character(sample))
  cf <- 1 - (NMLG/sample)
  cf <- (sample/(sample - 1))*cf
  return(cf)
}


datalist <- full_join(datalist, divlist, by = "pop") %>%
  full_join(loclist, by = "pop") %>%
  mutate(CF = CF(NMLG, sample))

vals <- datalist %>%
  mutate(mean.rd = vapply(samples.rd, mean, numeric(1), na.rm = TRUE)) %>%
  mutate(sd.rd = vapply(samples.rd, sd, numeric(1), na.rm = TRUE)) %>%
  mutate(var.rd = vapply(samples.rd, var, numeric(1), na.rm = TRUE)) %>%
  mutate(min.rd = vapply(samples.rd, min, numeric(1), na.rm = TRUE)) %>%
  mutate(max.rd = vapply(samples.rd, max, numeric(1), na.rm = TRUE)) %>%
  select(-starts_with("samples"))

ggplot(vals, aes(x = sexrate, y = rbarD, fill = sample)) +
  geom_boxplot()

ggplot(vals, aes(x = Hexp, y = rbarD, color = log(p.rD))) +
  geom_point(alpha = 0.25) +
  scale_color_viridis()



ggplot(vals, aes(x = Hexp, y = rbarD, color = log(p.rD))) +
  geom_point(alpha = 0.25) +
  scale_color_viridis() +
  facet_grid(sexrate~sample)


ggplot(vals, aes(x = Evenness, y = E.5, color = rbarD)) +
  geom_point(alpha = 0.25) +
  scale_color_viridis() +
  facet_grid(sexrate~sample)


ggplot(filter(vals, p.rD > 0.75), aes(x = Evenness, y = rbarD, color = E.5)) +
  geom_point(alpha = 0.25) +
  scale_color_viridis(option = "magma", direction = -1) +
  facet_grid(sexrate~sample) +
  theme_dark()


ggplot(vals, aes(x = sexrate, y = E.5, fill = sample)) +
  geom_boxplot()
ggplot(vals, aes(x = sexrate, y = H, fill = sample)) +
  geom_boxplot()
ggplot(vals, aes(x = sexrate, y = B, fill = sample)) +
  geom_boxplot()

ggplot(vals, aes(x = sexrate, y = E.5, color = log(p.rD))) +
  geom_jitter(alpha = 0.25, height = 0) +
  geom_boxplot(color = "black", width = 0.5, alpha = 0.25) +
  scale_color_viridis(option = "viridis",
                      breaks = log(c(0.005, 0.01, 0.025, 0.05, 0.1)),
                      labels = c(0.005, 0.01, 0.025, 0.05, 0.1)
                      ) +
  facet_wrap(~sample)

ggplot(vals, aes(x = CF, y = B, color = sample)) +
  geom_point() +
  facet_wrap(~sample)

ggplot(vals, aes(x = uSimp.var, y = uSimp)) +
  geom_point(aes(color = uSimp.var), position = "jitter", alpha = 0.25) +
  scale_color_viridis(option = "viridis") +
  facet_wrap(~sample)

ggplot(vals, aes(x = sexrate, y = p.rD, fill = sample)) +
  geom_boxplot() +
  scale_y_log10(breaks = c(0.01, 0.1, 1),
                minor_breaks = c((1:9)*0.001, (2:9)*0.01, (2:9)*0.1))

ggplot(vals, aes(x = sexrate, y = rbarD, color = log(p.rD))) +
  geom_jitter(alpha = 0.25, height = 0) +
  geom_boxplot(color = "black", notch = TRUE, width = 0.5, alpha = 0.25) +
  scale_color_viridis(option = "viridis", breaks = log(c(0.005, 0.01, 0.025, 0.05, 0.1)), labels = c(0.005, 0.01, 0.025, 0.05, 0.1)) +
  facet_wrap(~sample, nrow = 1)

ggplot(vals, aes(x = sexrate, y = rbarD - rbarDcc, color = log(p.rD))) +
  geom_point(alpha = 0.25, position = "jitter") +
  geom_boxplot(color = "black", notch = TRUE, width = 0.5, alpha = 0.25) +
  scale_color_viridis(option = "viridis", breaks = log(c(0.005, 0.01, 0.025, 0.05, 0.1)), labels = c(0.005, 0.01, 0.025, 0.05, 0.1)) +
  facet_wrap(~sample, nrow = 1)

ggplot(vals, aes(x = rbarD, y = rbarDcc, color = log(p.rD))) +
  geom_point(alpha = 0.25) +
  scale_color_viridis(option = "viridis", breaks = log(c(0.005, 0.01, 0.025, 0.05, 0.1)), labels = c(0.005, 0.01, 0.025, 0.05, 0.1)) +
  geom_abline(slope = 1) +
  facet_grid(sample~sexrate) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  xlab(expression(bar(r)[d])) +
  ylab(expression(paste(bar(r)[d], " (clone-corrected)")))

