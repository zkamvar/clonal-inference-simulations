sims <- getfile()
sim <- read.table(sims$files, header = TRUE, sep = " ")
sim$Generation <- as.numeric(vapply(strsplit(as.character(sim$File), "_"), function(x) strsplit(x[10], "\\.")[[1]][1], "butt"))
sim$Samp.Size <- factor(sim$Samp.Size, levels = c(10,25,50,100))

simnarm <- sim[-which(is.na(sim$rbarD)), ]

Ia.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation*1000), y = Ia, fill = Samp.Size)) + facet_wrap(~Sex.Rate) + theme_bw() + scale_fill_grey(start = 0.4, end = 1)
rbarD.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation*1000), y = rbarD, fill = Samp.Size)) + facet_wrap(~Sex.Rate) + theme_bw() + scale_fill_grey(start = 0.4, end = 1)


Ia.point <- ggplot(simnarm, aes(x = factor(Generation*1000), y = Ia, color = Samp.Size, group = Samp.Size, alpha = 0.5)) + geom_point(position = "jitter") + stat_smooth() + facet_grid(Sex.Rate ~ Samp.Size) + theme_bw()
rbarD.point <- ggplot(simnarm, aes(x = factor(Generation*1000), y = rbarD, color = Samp.Size, group = Samp.Size, alpha = 0.5)) + geom_point(position = "jitter") + stat_smooth() + facet_grid(Sex.Rate ~ Samp.Size) + theme_bw()
